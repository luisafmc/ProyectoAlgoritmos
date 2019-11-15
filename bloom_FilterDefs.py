import pyhash
import sys
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt1

bloomFilterSize=10
bit_vector=[]


#hashFunctions
fnv=pyhash.fnv1a_32()
mur=pyhash.murmur3_32()
lookup=pyhash.lookup3()
super1=pyhash.super_fast_hash()
city=pyhash.city_64()
spooky=pyhash.spooky_32()
farm=pyhash.farm_32()
metro=pyhash.metro_64()
mum=pyhash.mum_64()
xx=pyhash.xx_32()

#10 hash functions
hashfuncs=[fnv,mur,lookup,super1,city,spooky,farm,metro,mum,xx]

#Create the bloom filter and add the kmers to bit vector
def insertBloom(kmer, hashFuncCount):
    global bloomFilterSize
    global bit_vector
    index=0
    for hf in hashfuncs:
        if(index<=hashFuncCount):
            if(bit_vector[hf(kmer)%bloomFilterSize]==0):
                for hf2 in hashfuncs:
                    bit_vector[hf2(kmer)%bloomFilterSize]=1
            index+=1
        else:
            index+=1
            break
        #print("not avaliable")

#Inspect the filter
def lookFilter(kmer,hashFuncCount):
    global bit_vector
    global bloomFilterSize
    index=0
    for hs in hashfuncs:
        if(bit_vector[hs(kmer)%bloomFilterSize]==1):
            return True
        if(index>=hashFuncCount):
            break
        index+=1
    return False
            
            
#Calculate the reverse complement
def reverse(seq):
    """Returns a reversed string"""
    return seq[::-1]


def complement(seq):
    """Returns a complement DNA sequence"""
    complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    seq_list = list(seq)
    seq_list = [complement_dict[base] for base in seq_list]
    return ''.join(seq_list)


def reverse_complement(seq):
    """"Returns a reverse complement DNA sequence"""
    seq = reverse(seq)
    seq = complement(seq)
    return seq

#Principal function
def main(argv):
    #Bloom filter size
    global bloomFilterSize
    global bit_vector
    bloomFilterSize=int(argv[0])
    
    #BloomFilter
    bit_vector=[0]*bloomFilterSize
        
    #HashTable
    hashT={}
    
    #KmerSize
    k=int(argv[1])
    
    #filename
    filename=argv[2]
    string = ''

    #HashFunctions to use
    hashFuncC=int(argv[3])

    #File output
    name=filename.split("_")[-2]

    fileOutput = open(str(name)+"_"+str(k)+"_out.fasta" , "w")

    print("bit_vector_size",len(bit_vector))
    print("Procesing " + filename)
    print("The number of hash functions in use are " + str(hashFuncC) + "\n")

    f1=open(filename,'r')
    try:
        #print("lectura de archiv")
        index = 0        
        for l1 in f1:
            if l1.startswith('>'):
                if (index == 0):
                    index +=1
                else:
                    string+=l1.rstrip("\r\n")
            else:
                if (index == 0):
                    index +=1
                else:
                    string+=l1.rstrip("\r\n")

    finally:
        f1.close()

  


    kmers=[]
    idx=0
    while(len(string)-idx>=k):
        kmers.append(string[idx:idx+k].rstrip("\r\n"))
        idx+=1

    #print("---kmers---")
    #print(len(kmers))

    for kmer in kmers:
        xrep=''
        x=kmer
        xrep=x
        #xrev=reverse_complement(kmer)
        #if(xrev>x):
        #    xrep=x
        #else:
        #    xrep=xrev
        if(lookFilter(xrep,hashFuncC)):
            if(xrep not in hashT.keys()):
                hashT[xrep]=0
        else:
            insertBloom(xrep,hashFuncC)
    for kmer in kmers:
        xrep=''
        x=kmer
        xrep=x
        #xrev=reverse_complement(kmer)
        #if(xrev>x):
        #    xrep=x
        #else:
        #    xrep=xrev
        if(xrep in hashT.keys()):
            hashT[xrep]+=1

    definiteDict={}
    uniqueUnfiltered=0

    for key in hashT.keys():
        if(hashT[key]!=1):
            definiteDict[key]=hashT[key]
        else:
            uniqueUnfiltered+=1
            
    #print (definiteDict)    
    #print(uniqueUnfiltered)
    Max=max(definiteDict.values())
    Len=len(definiteDict.keys())

	#Plot of frequency (normal)
    plt1.hist(hashT.values(),bins=100,alpha=1,color="orange")
    plt1.xlabel('Coverage')
    plt1.ylabel('Frequency')
    plt1.title(r'Normal Histogram of kmers profile')
    plt1.savefig(str(name) + "_" + str(k)+ "_" + str(bloomFilterSize)+ "_" +"Normal_hist.png",dpi=200)
    plt1.close()
    
    #Plot of frequency (bloom filter)
    plt.hist(definiteDict.values(),bins=100,alpha=1,color="blue")
    plt.xlabel('Coverage')
    plt.ylabel('Frequency')
    plt.title(r'Bloom filter Histogram of kmers profile')
    plt.savefig(str(name) + "_" + str(k)+ "_" + str(bloomFilterSize) + "_" +"Bloomfilter_hist.png",dpi=200)
    plt.close()

    #Fasta generator
    for k2, v in definiteDict.items():
        ###Output the header
        fileOutput.write(">" + str(v) + "\n")
        fileOutput.write(str(k2) + "\n")
    
    fileOutput.close()
    print ("Fasta file done. Name is " + str(str(name)+"_"+str(k)+"_out.fasta") + "\n")

    #Statistics

    
    print ("Unique: " + str(uniqueUnfiltered))
    print ("Distinct: " + str(len(hashT.keys())))
    print ("Total: " + str(len(kmers)))
    print ("Max_count: " + str(Max) + "\n")


# main function
if __name__ == "__main__":
   main(sys.argv[1:])
    
