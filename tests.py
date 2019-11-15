#tests
import subprocess
import os
cmd='cmd C'

b0=5000
bf=100000
saltos=5000

hfo=1
sep=";"
hff=5

print(os.getcwd())
f=open("results_15_fly.csv","w")
f.write("HashFunctions;BloomFilterSize;UniqueKmers;DistinctKmers")
try:
	while(b0<bf):
		while(hfo<hff):
			data, temp=os.pipe()

			os.write(temp,bytes("5 10\n","utf-8"))
			os.close(temp)

			s=subprocess.check_output("python bloom_FilterDefs.py "+str(b0)+ " 15 fasta/GCF_000001215.4_D.melanogaster_genomic.fna "+str(hff),shell=True)	
			#print(s.decode("utf-8"))

			op=s.decode("utf-8").split("\n")
			#print(op[2])

			hfsT=op[2].split("The number of hash functions in use are ")[1]
			#print(op[0])
			bfs=op[0].split("bit_vector_size ")[1]
			#print(op[6])
			ukc=op[6].split("Unique: ")[1]
			#print(op[7])
			dkc=op[7].split("Distinct: ")[1]

			f.write(hfsT+sep+bfs+sep+ukc+sep+dkc+"\n")
			print()
			hfo+=1
		b0+=saltos
		hfo=1
		print(b0)
finally:
	f.close()



