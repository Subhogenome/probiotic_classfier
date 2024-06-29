import argparse
from Bio import SeqIO
import numpy as np
import pickle
from io import StringIO
from itertools import product
import sys
import warnings
import collections
import lightgbm as lgb
import os.path
import pandas as pd
import streamlit as st
import sklearn
magicEnabled = True
loaded_model = pickle.load(open('finalized_model.sav', 'rb'))
if os.path.exists("dataset.csv"):
    os.remove("dataset.csv")


if os.path.exists("output.fasta"):
    os.remove("output.fasta")
if os.path.exists("Results.csv"):
    os.remove("Results.csv")
    
if os.path.exists("myfile.txt"):
    os.remove("myfile.txt")
    

open('myfile.txt', 'w')   

    
stro="output.fasta"
finput=open(stro,'w+') 
foutput=str("dataset.csv")

labelDataset = str("DNA")
ksize = int(8)
seq = int(1)
stepw = 1

foutput=str("dataset.csv")
poutput=str("output.fasta")

if seq == 1:
    	caracteres = ["A", "C", "G", "T"]
else:
    	caracteres = ["A", "C", "G", "U"]


  
def preprocessing(pinput,poutput):
    alphabet = ("B|D|E|F|H|I|J|K|L|M|N|O|P|Q|R|S|V|W|X|Y|Z")
    file = open(poutput, 'a')
    for seq_record in SeqIO.parse(pinput, "fasta"):
        name_seq = seq_record.name
        seq = seq_record.seq
        des=seq_record.description
        des= des.split()
        speices=des[1]
        name = des[2]
        Genome=speices+" "+name
        file1 = open("myfile.txt","w")
        file1.write(Genome)
        file1.close() 
        if re.search(alphabet, str(seq)) is not None:
            print(name_seq)
            print("Removed Sequence")
        else:
            file.write(">%s" % (str(name_seq)))
            file.write("\n")
            file.write(str(seq))
            file.write("\n")
            print(name_seq)
            print("Included Sequence")
    print("Finished")


def preprocessingD(pinput,poutput):
    alphabet = ("B|D|E|F|H|I|J|K|L|M|N|O|P|Q|R|S|V|W|X|Y|Z")
    file = open(poutput, 'a')
    for seq_record in SeqIO.parse(pinput, "fasta"):
        name_seq = seq_record.name
        seq = seq_record.seq
        

        if re.search(alphabet, str(seq)) is not None:
            st.write(name_seq)
            st.write("Removed Sequence")
        else:
            file.write(">%s" % (str(name_seq)))
            file.write("\n")
            file.write(str(seq))
            file.write("\n")
            print(name_seq)
            print("Included Sequence")
    print("Finished")
    pinput.close()
import sys
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import collections
import argparse
from Bio import SeqIO
from itertools import product


def perms():
    global listTotal, caracteres
    listTotal = []
    file = open(foutput, 'a')
    file.write("%s," % ("nameseq"))
    for k in range(1, ksize+1):
        permsList = [''.join(str(i) for i in x) for x in product(caracteres, repeat=k)]
        # print (permsList)
        for perm in permsList:
            # print (perm)
            listTotal.append(perm)
            file.write("%s," % (str(perm)))
    file.write("label")
    file.write("\n")
    return


def chunks(seq, win, step):
    seqlen = len(seq)
    for i in range(0,seqlen,step):
        j = seqlen if i+win>seqlen else i+win
        yield seq[i:j]
        if j==seqlen: break
    return        
    

def chunksTwo(seq, win):
    seqlen = len(seq)
    for i in range(seqlen):
        j = seqlen if i+win>seqlen else i+win
        yield seq[i:j]
        if j==seqlen: break
    return

            
def fileRecord(name_seq):
    file = open(foutput, 'a')
    file.write("%s," % (name_seq))
    for x in probabilities:
        # print (x)
        file.write("%s," % (str(x[1])))
    file.write(labelDataset)
    file.write("\n")
    print("Recorded Sequence: %s" % (name_seq))
    return
    

def findKmers():
    global probabilities
    perms()
    for seq_record in SeqIO.parse(finput, "fasta"):
        seq = seq_record.seq
        seq = seq.upper()	
        name_seq = seq_record.name
        probabilities = []
        for k in range(1, ksize+1):
            kmer = {}
            totalWindows = (len(seq) - k) + 1 # (L - k + 1)
            permsList = [''.join(str(i) for i in x) for x in product(caracteres, repeat=k)]
            for key in permsList:
                kmer[key] = 0
            kmer = collections.OrderedDict(sorted(kmer.items()))
            for subseq in chunksTwo(seq, k):
                if subseq in kmer:
                    # print(subseq)
                    kmer[subseq] = kmer[subseq] + 1
                else:
                    kmer[subseq] = 1
            for key, value in kmer.items():
                # print (key)
                # print (value)
                probabilities.append([str(key), value/totalWindows])
        fileRecord(name_seq)
    return

def convert_df(df):
   return df.to_csv().encode('utf-8')

def prediction(input_):
   file1 = open("myfile.txt","r+")
   Genome=file1.read()
   input_data=pd.read_csv("dataset.csv")
   X=input_data.iloc[:,2:87380]
   prediction = pd.DataFrame(loaded_model.predict_proba(X))
   prediction.rename(columns = {0:'Probiotic', 1:'Pathogenic'}, inplace = True)
   input_data['Pathogenic'] = pd.DataFrame(prediction, columns=['Pathogenic'])
   input_data['Probiotic'] = pd.DataFrame(prediction, columns=['Probiotic'])
   results=input_data.iloc[:,[0,87382,87383]]
   results=results.set_index('nameseq')
   results.loc[Genome+" "+"Prediction"] = results.mean()
   data={'Name':[Genome, Genome],"Probablites":[results["Pathogenic"].mean(),results["Probiotic"].mean()],"Type":["Non-Probiotic",'Probiotic']}
   df=pd.DataFrame(data)
   fig = px.bar(df, x="Probablites", y="Name", color='Type', orientation='h',
             hover_data=["Name"],
             height=400,title=(Genome+" "+"Overall Prediction"))
   with st.container():
     st.plotly_chart(fig)
     st.dataframe(results)
     st.download_button("Press to Download",convert_df(results),"Result.csv","text/csv",key='download-csv')

def predictionD(input_):

   input_data=pd.read_csv("dataset.csv",error_bad_lines=False)
   X=input_data.iloc[:,2:87380]
   prediction = pd.DataFrame(loaded_model.predict_proba(X))
   prediction.rename(columns = {0:'Probiotic', 1:'Pathogenic'}, inplace = True)
   input_data['Pathogenic'] = pd.DataFrame(prediction, columns=['Pathogenic'])
   input_data['Probiotic'] = pd.DataFrame(prediction, columns=['Probiotic'])
   results=input_data.iloc[:,[0,87382,87383]]
   results=results.set_index('nameseq')
   results.loc["Overall Prediction"] = results.mean()
   data={'Name':["Overall Prediction" , "Overall Prediction"],"Probablites":[results["Pathogenic"].mean(),results["Probiotic"].mean()],"Type":["Non-Probiotic",'Probiotic']}
   df=pd.DataFrame(data)
   fig = px.bar(df, x="Probablites", y="Name", color='Type', orientation='h',
             hover_data=["Name"],
             height=400,title=("Overall Prediction"))
   with st.container():
     st.plotly_chart(fig)
     st.dataframe(results)
     st.download_button("Press to Download",convert_df(results),"Result.csv","text/csv",key='download-csv')


def argmine():
    cmd = 'amrfinder --nucleotide output.fasta --output output.txt'
    os.system(cmd)
    data=pd.read_csv('output.txt',sep='\t')
    st.dataframe(data)







def main():

    st.image("logo.png")


    
     # giving a title
    st.title('e-predict')
    finput = st.file_uploader("upload file", type={".fasta","fna"})
    genre = st.radio(
    "Specify the type of genome being used",
    ('Draft Genome', 'Complete sequence'))

    if finput is not None:

            # To convert to a string based IO:
            stringio = StringIO(finput.getvalue().decode("utf-8"))

    
    
    # code for Prediction
    diagnosis = ''
    
    
    
    if st.button('Predict results'):
        if genre == 'Draft Genome':
           
            st.info("loading Genome")
            preprocessingD(stringio,poutput)
            st.info("Creating K-mers")
        
            findKmers()
            st.info("Predicting Results")
            diagnosis = predictionD(stringio)
            argmine()
        
        else:
          st.info("loading Genome")
          preprocessing(stringio,poutput)
          st.info("Creating K-mers")

        
          findKmers()
          st.info("Predicting Results")
          diagnosis = prediction(stringio)
          st.info("AMRfinder results")
          #argmine()


if __name__ == '__main__':
    main()                           






