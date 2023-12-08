"""
.. module:: MRKmeansDef

MRKmeansDef
*************

:Description: MRKmeansDef

    

:Authors: bejar
    

:Version: 

:Created on: 17/07/2017 7:42 

"""

from functools import reduce
from mrjob.job import MRJob
from mrjob.step import MRStep
import collections
import itertools

__author__ = 'bejar'


class MRKmeansStep(MRJob):
    prototypes = {}

    def jaccard(self, prot, doc):
        """
        Compute here the Jaccard similarity between  a prototype and a document
        prot should be a list of pairs (word, probability)
        doc should be a list of words
        Words must be alphabeticaly ordered

        The result should be always a value in the range [0,1]
        """
        w1 = 0   #prot
        w2 = 0   #doc
        dot = 0

        pointerw1 = 0
        pointerw2 = 0

        d1 = False
        d2 = False


        while not (d1 and d2):

            if (prot[pointerw1][0] < doc[pointerw2] and not d1) or(not d1 and d2):   #pointer of the first word is on a smaller word
                w1 += pow(prot[pointerw1][1],2)  #inner sum
                pointerw1 += 1  #move pointer
        
            elif (prot[pointerw1][0] > doc[pointerw2] and not d2) or(not d2 and d1): #pointer of the second word is on a smaller word
                w2 += 1
                pointerw2 += 1
    
            else: #both pointers are on the same word
                w1 += pow(prot[pointerw1][1],2) #inner sum
                w2 += 1
                dot += prot[pointerw1][1]  #append to dot product
                pointerw1 += 1 #move both pointers
                pointerw2 += 1

            if pointerw1 == len(prot):  #checks of all words in a vector are read and prevents that pointer goes to far
                pointerw1 -= 1
                d1 = True

            if pointerw2 == len(doc):
                pointerw2 -= 1
                d2 = True
            
        return dot/(w1+w2-dot) 


    def configure_args(self):
        """
        Additional configuration flag to get the prototypes files

        :return:
        """
        super(MRKmeansStep, self).configure_args()
        self.add_file_arg('--prot')

    def load_data(self):
        """
        Loads the current cluster prototypes

        :return:
        """
        f = open(self.options.prot, 'r')
        for line in f:
            cluster, words = line.split(':')
            cp = []
            for word in words.split():
                cp.append((word.split('+')[0], float(word.split('+')[1])))
            self.prototypes[cluster] = cp

    def assign_prototype(self, _, line):
        """
        This is the mapper it should compute the closest prototype to a document

        Words should be sorted alphabetically in the prototypes and the documents

        This function has to return at list of pairs (prototype_id, document words)

        You can add also more elements to the value element, for example the document_id
        """

        # Each line is a string docid:wor1 word2 ... wordn
        doc, words = line.split(':')
        lwords = words.split()

        min_prot = list(self.prototypes.keys())[0]
        val_prot = 0

        for p in self.prototypes:

            cur_val = self.jaccard(self.prototypes[p], lwords)
            if cur_val > val_prot:
                min_prot=p
                val_prot = cur_val

        # Return pair key, value
        yield min_prot, line         #key: prot_id (cluster), value: docid:wor1 word2 ... wordn

    def aggregate_prototype(self, key, values):
        """
        input is cluster and all the documents it has assigned
        Outputs should be at least a pair (cluster, new prototype)

        It should receive a list with all the words of the documents assigned for a cluster

        The value for each word has to be the frequency of the word divided by the number
        of documents assigned to the cluster

        Words are ordered alphabetically but you will have to use an efficient structure to
        compute the frequency of each word

        :param key:
        :param values:
        :return:
        """

        doc_num = len(values)

        word_list = []
        doc_list = []

        for i in values:
            doc, words = i.split(':')
            word_list += words.split()
            doc_list.append(doc)


        counter = dict(collections.Counter(word_list))

        res = []

        for i in counter:
            res.append((i, counter[key]/doc_num))
            
        yield key, (res, doc_list)

    def steps(self):
        return [MRStep(mapper_init=self.load_data, mapper=self.assign_prototype,
                       reducer=self.aggregate_prototype)
            ]


if __name__ == '__main__':
    MRKmeansStep.run()