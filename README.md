# AAV integration detection script

Inherited from Laurence Wilson.  He wrote it to detect AAV integration sites in data from mice (FRG or wild-type) treated with AAV, from collaborators at Westmead (CMRI).

## 20-05-2019

Claus at CMRI suggested possible change: he noticed that some of the human part of some of  the reads aligns equally well to multiple parts of the human genome.  These are likely parts of the genome where duplications have occurred  He argues that an integration site can't be uniquely defined if the reads can't be localised to a unique location in the genome.  He wanted some way of flagging up these reads.

