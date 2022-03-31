

msa = "rep_and_twist.msa"
name_to_seq = {name:seq for name,seq in sfp(open(msa))}
refseq = name_to_seq["MN908947.3"]


index_to_pos = {}
pos_to_indexes = defaultdict(list)

# the add1 is necessary for consistency with variants build under:
# when looking at variant we concat all base observed until the pos for which ref is != "-". 
# we then use this pos for the concat bases.
# here we use the variable add1 for the same purpose
pos=0
add1=0
for index,nuc in enumerate(refseq):
	if (nuc=="-")&(add1==0):
		add1 = 1
	if nuc!="-":
		pos+=1
		add1=0
	index_to_pos[index]=pos+add1
	pos_to_indexes[pos+add1].append(index)


trim = [300,29700]
Variants = []
for name,seq in name_to_seq.items():
	temp=""
	for index,nuc in enumerate(refseq):
		pos = index_to_pos[index]
		if trim[0]<=pos<=trim[1]:
			var_nuc = seq[index]
			if temp:
				temp+=var_nuc
				if nuc!="-":
					if var_nuc==nuc:
						Variants.append([name,index_to_pos[index],temp])
						temp=""
			else:
				if nuc!="-":
					if nuc!=var_nuc:
						Variants.append([name,index_to_pos[index],var_nuc])
				else:
					if nuc!=var_nuc:
						temp+=var_nuc


with open("db_snv.tsv","w") as handle:
	handle.writelines("%s\t%s\t%s\n"%(name,str(pos),nuc) for name,pos,nuc in Variants)

