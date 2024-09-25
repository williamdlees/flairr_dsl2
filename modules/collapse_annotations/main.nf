

process collapse_annotations {

	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${outfile}+passed.tsv$/) "collapsed/${name}_collapsed_passed_${alignment_suffix}.tsv"}
	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${outfile}+failed.*$/) "collapsed/${name}_collapsed_failed_${alignment_suffix}.tsv"}
	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename == ".command.out") "collapsed/${name}_collapsed_${alignment_suffix}.log"}
	
	input:
		path airrFile
		val alignment_suffix

	output:
		path("${outfile}" + "passed.tsv"), emit: output
		path("${outfile}" + "failed*") optional true
		path(".command.out") optional true

	script:
		name = params.sample_name
		conscount_min = params.collapse_annotations.conscount_min
		n_max = params.collapse_annotations.n_max


		outfile = airrFile.toString() - '.tsv' + alignment_suffix + "_collapsed-"

		if(airrFile.getName().endsWith(".tsv")) {	
			"""
			#!/usr/bin/env python3
			
			from pprint import pprint
			from collections import OrderedDict,Counter
			import itertools as it
			import datetime
			import pandas as pd
			import glob, os
			import numpy as np
			import re
			
			# column types default
			
			# dtype_default={'junction_length': 'Int64', 'np1_length': 'Int64', 'np2_length': 'Int64', 'v_sequence_start': 'Int64', 'v_sequence_end': 'Int64', 'v_germline_start': 'Int64', 'v_germline_end': 'Int64', 'd_sequence_start': 'Int64', 'd_sequence_end': 'Int64', 'd_germline_start': 'Int64', 'd_germline_end': 'Int64', 'j_sequence_start': 'Int64', 'j_sequence_end': 'Int64', 'j_germline_start': 'Int64', 'j_germline_end': 'Int64', 'v_score': 'Int64', 'v_identity': 'Int64', 'v_support': 'Int64', 'd_score': 'Int64', 'd_identity': 'Int64', 'd_support': 'Int64', 'j_score': 'Int64', 'j_identity': 'Int64', 'j_support': 'Int64'}
			
			SPLITSIZE=2
			
			class CollapseDict:
				def __init__(self,iterable=(),_depth=0,
							 nlim=10,conscount_flag=False):
					self.lowqual={}
					self.seqs = {}
					self.children = {}
					self.depth=_depth
					self.nlim=nlim
					self.conscount=conscount_flag
					for fseq in iterable:
						self.add(fseq)
			
				def split(self):
					newseqs = {}
					for seq in self.seqs:
						if len(seq)==self.depth:
							newseqs[seq]=self.seqs[seq]
						else:
							if seq[self.depth] not in self.children:
								self.children[seq[self.depth]] = CollapseDict(_depth=self.depth+1)
							self.children[seq[self.depth]].add(self.seqs[seq],seq)
					self.seqs=newseqs
			
				def add(self,fseq,key=None):
					#if 'duplicate_count' not in fseq: fseq['duplicate_count']='1'
					if 'KEY' not in fseq:
						fseq['KEY']=fseq['sequence_vdj'].replace('-','').replace('.','')
					if 'ISOTYPECOUNTER' not in fseq:
						fseq['ISOTYPECOUNTER']=Counter([fseq['c_call']])
					if 'VGENECOUNTER' not in fseq:
						fseq['VGENECOUNTER']=Counter([fseq['v_call']])
					if 'JGENECOUNTER' not in fseq:
						fseq['JGENECOUNTER']=Counter([fseq['j_call']])
					if key is None:
						key=fseq['KEY']
					if self.depth==0:
						if (not fseq['j_call'] or not fseq['v_call']):
							return
						if fseq['sequence_vdj'].count('N')>self.nlim:
							if key in self.lowqual:
								self.lowqual[key] = combine(self.lowqual[key],fseq,self.conscount)
							else:
								self.lowqual[key] = fseq
							return
					if len(self.seqs)>SPLITSIZE:
						self.split()
					if key in self.seqs:
						self.seqs[key] = combine(self.seqs[key],fseq,self.conscount)
					elif (self.children is not None and
						  len(key)>self.depth and
						  key[self.depth] in self.children):
						self.children[key[self.depth]].add(fseq,key)
					else:
						self.seqs[key] = fseq
			
				def __iter__(self):
					yield from self.seqs.items()
					for d in self.children.values():
						yield from d
					yield from self.lowqual.items()
			
				def neighbors(self,seq):
					def nfil(x): return similar(seq,x)
					yield from filter(nfil,self.seqs)
					if len(seq)>self.depth:
						for d in [self.children[c]
								  for c in self.children
								  if c=='N' or seq[self.depth]=='N' or c==seq[self.depth]]:
							yield from d.neighbors(seq)
			
				def fixedseqs(self):
					return self
					ncd = CollapseDict()
					for seq,fseq in self:
						newseq=seq
						if 'N' in seq:
							newseq=consensus(seq,self.neighbors(seq))
							fseq['KEY']=newseq
						ncd.add(fseq,newseq)
					return ncd
			
			
				def __len__(self):
					return len(self.seqs)+sum(len(c) for c in self.children.values())+len(self.lowqual)
			
			def combine(f1,f2, conscount_flag):
				def val(f): return -f['KEY'].count('N'),(int(f['consensus_count']) if 'consensus_count' in f else 0)
				targ = (f1 if val(f1) >= val(f2) else f2).copy()
				if conscount_flag:
					targ['consensus_count'] =  int(f1['consensus_count'])+int(f2['consensus_count'])
				targ['duplicate_count'] =  int(f1['duplicate_count'])+int(f2['duplicate_count'])
				targ['ISOTYPECOUNTER'] = f1['ISOTYPECOUNTER']+f2['ISOTYPECOUNTER']
				targ['VGENECOUNTER'] = f1['VGENECOUNTER']+f2['VGENECOUNTER']
				targ['JGENECOUNTER'] = f1['JGENECOUNTER']+f2['JGENECOUNTER']
				return targ
			
			def similar(s1,s2):
				return len(s1)==len(s2) and all((n1==n2 or n1=='N' or n2=='N')
											  for n1,n2 in zip(s1,s2))
			
			def basecon(bases):
				bases.discard('N')
				if len(bases)==1: return bases.pop()
				else: return 'N'
			
			def consensus(seq,A):
				return ''.join((basecon(set(B)) if s=='N' else s) for (s,B) in zip(seq,zip(*A)))
				
			def unpack_counter(c):
				return ', '.join([f'{x[0]}: {x[1]}' for x in [(k, c[k]) for k in sorted(c.keys())]])
			
			n_lim = int('${n_max}')
			conscount_filter = int('${conscount_min}')
			
			df = pd.read_csv('${airrFile}', sep = '\t') #, dtype = dtype_default)
			
			# make sure that all columns are int64 for createGermline
			idx_col = df.columns.get_loc("cdr3")
			cols =  [col for col in df.iloc[:,0:idx_col].select_dtypes('float64').columns.values.tolist() if not re.search('support|score|identity|freq', col)]
			df[cols] = df[cols].apply(lambda x: pd.to_numeric(x.replace("<NA>",np.NaN), errors = "coerce").astype("Int64"))
			
			conscount_flag = False
			if 'consensus_count' in df: conscount_flag = True
			if not 'duplicate_count' in df:
				df['duplicate_count'] = 1
			if not 'c_call' in df or not 'isotype' in df or not 'prcons' in df or not 'primer' in df or not 'reverse_primer' in df:
				if 'c_call' in df:
					df['c_call'] = df['c_call']
				elif 'isotype' in df:
					df['c_call'] = df['isotype']
				elif 'primer' in df:
					df['c_call'] = df['primer']
				elif 'reverse_primer' in df:
					df['c_call'] = df['reverse_primer']    
				elif 'prcons' in df:
					df['c_call'] = df['prcons']
				elif 'barcode' in df:
					df['c_call'] = df['barcode']
				else:
					df['c_call'] = 'Ig'
			
			# removing sequenes with duplicated sequence id    
			dup_n = df[df.columns[0]].count()
			df = df.drop_duplicates(subset='sequence_id', keep='first')
			dup_n = str(dup_n - df[df.columns[0]].count())
			df['c_call'] = df['c_call'].astype('str').replace('<NA>','Ig')
			#df['consensus_count'].fillna(2, inplace=True)
			nrow_i = df[df.columns[0]].count()
			df = df[df.apply(lambda x: x['sequence_alignment'][0:(x['v_germline_end']-1)].count('N')<=n_lim, axis = 1)]
			low_n = str(nrow_i-df[df.columns[0]].count())
			
			df['sequence_vdj'] = df.apply(lambda x: x['sequence_alignment'].replace('-','').replace('.',''), axis = 1)
			header=list(df.columns)
			fasta_ = df.to_dict(orient='records')
			c = CollapseDict(fasta_,nlim=10)
			d=c.fixedseqs()
			header.append('ISOTYPECOUNTER')
			header.append('VGENECOUNTER')
			header.append('JGENECOUNTER')
			
			rec_list = []
			for i, f in enumerate(d):
				rec=f[1]
				rec['sequence']=rec['KEY']
				rec['consensus_count']=int(rec['consensus_count']) if conscount_flag else None
				rec['duplicate_count']=int(rec['duplicate_count'])
				rec['ISOTYPECOUNTER'] = unpack_counter(rec['ISOTYPECOUNTER'])
				rec['VGENECOUNTER'] = unpack_counter(rec['VGENECOUNTER'])
				rec['JGENECOUNTER'] = unpack_counter(rec['JGENECOUNTER'])
				rec_list.append(rec)
			df2 = pd.DataFrame(rec_list, columns = header)        
			
			df2 = df2.drop('sequence_vdj', axis=1)
			
			collapse_n = str(df[df.columns[0]].count()-df2[df2.columns[0]].count())

			# removing sequences without J assignment and non functional
			nrow_i = df2[df2.columns[0]].count()
			cond = (~df2['j_call'].str.contains('J')|df2['productive'].isin(['F','FALSE','False']))
			df_non = df2[cond]
			
			
			df2 = df2[df2['productive'].isin(['T','TRUE','True'])]
			cond = ~(df2['j_call'].str.contains('J'))
			df2 = df2.drop(df2[cond].index.values)
			
			non_n = nrow_i-df2[df2.columns[0]].count()
			#if conscount_flag:
			#   df2['consensus_count'] = df2['consensus_count'].replace(1,2)
			
			# removing sequences with low cons count
			
			filter_column = "duplicate_count"
			if conscount_flag: filter_column = "consensus_count"
			df_cons_low = df2[df2[filter_column]<conscount_filter]
			nrow_i = df2[df2.columns[0]].count()
			df2 = df2[df2[filter_column]>=conscount_filter]
			
			
			cons_n = str(nrow_i-df2[df2.columns[0]].count())
			nrow_i = df2[df2.columns[0]].count()    
			
			df2.to_csv('${outfile}'+'passed.tsv', sep = '\t',index=False) #, compression='gzip'
			
			pd.concat([df_cons_low,df_non]).to_csv('${outfile}'+'failed.tsv', sep = '\t',index=False)
			
			print(str(low_n)+' Sequences had N count over 10')
			print(str(dup_n)+' Sequences had a duplicated sequence id')
			print(str(collapse_n)+' Sequences were collapsed')
			print(str(df_non[df_non.columns[0]].count())+' Sequences were declared non functional or lacked a J assignment')
			#print(str(df_cons_low[df_cons_low.columns[0]].count())+' Sequences had a '+filter_column+' lower than threshold')
			print('Going forward with '+str(df2[df2.columns[0]].count())+' sequences')
			
			"""
		}else{
			"""
			"""
		}
}
