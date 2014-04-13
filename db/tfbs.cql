create (Fis:TransFactor:Protein {name:"Fis"}),(hns:TransFactor:Protein {name:"H-NS"}),
(fisBS1:TFBS:Feature {id:"ECK120011694",strand:"forward",start:3939411,end:3939425,pos:-121,seq:"tgtgattacgGTTGCAAAAACGGCAaattgcttgt",db:"RegulonDB"}),
(fisBS2:TFBS:Feature {id:"ECK120011692",strand:"forward",start:3939431,end:3939445,pos:-101,seq:"cggcaaattgCTTGTTTTATGGCACattaacgggg",db:"RegulonDB"}),
(fisBS3:TFBS:Feature {id:"ECK120011690",strand:"forward",start:3939462,end:3939476,pos:-70,seq:"ggggcttttgCTGAAAAAATGCGCGgtcagaaaat",db:"RegulonDB"}),
(fisBS4:TFBS:Feature {id:"ECK120011292",strand:"forward",start:4164240,end:4164254,pos:-143,seq:"gccaggagctGAACAATTATTGCCCgttttacagc",db:"RegulonDB"}),
(fisBS5:TFBS:Feature {id:"ECK120011290",strand:"forward",start:4164281,end:4164295,pos:-102,seq:"gcttcgaaacGCTCGAAAAACTGGCagttttaggc",db:"RegulonDB"}),
(fisBS6:TFBS:Feature {id:"ECK120011288",strand:"forward",start:4164312,end:4164326,pos:-71,seq:"aggctgatttGGTTGAATGTTGCGCggtcagaaaa",db:"RegulonDB"}),
(hnsBS1:TFBS:Feature {id:"ECK120011296",strand:"forward",start:4164300,end:4164314,pos:-83,seq:"actggcagttTTAGGCTGATTTGGTtgaatgttgc",db:"RegulonDB"}),
(hnsBS2:TFBS:Feature {id:"ECK120012498",strand:"forward",start:4164333,end:4164347,pos:-50,seq:"gcgcggtcagAAAATTATTTTAAATttcctcttgt",db:"RegulonDB"}),
(hnsBS3:TFBS:Feature {id:"ECK120012500",strand:"forward",start:4164357,end:4164371,pos:-26,seq:"tttcctcttgTCAGGCCGGAATAACtccctataat",db:"RegulonDB"})
with Fis,hns,fisBS1,fisBS2,fisBS3,fisBS4,fisBS5,fisBS6,hnsBS1,hnsBS2,hnsBS3
match (rrsBp1:Promoter),(rrsCp1:Promoter)
where rrsBp1.name="rrsBp1" and rrsCp1.name="rrsCp1"
create 
fisBS1-[:REGULATE {evidence:"Binding of purified proteins",evType:"strong",db:"RegulonDB",function:"activator"}]->rrsCp1,
fisBS2-[:REGULATE {evidence:"Binding of purified proteins",evType:"strong",db:"RegulonDB",function:"activator"}]->rrsCp1,
fisBS3-[:REGULATE {evidence:"Binding of purified proteins",evType:"strong",db:"RegulonDB",function:"activator"}]->rrsCp1,
fisBS4-[:REGULATE {evidence:"Binding of purified proteins",evType:"strong",db:"RegulonDB",function:"activator"}]->rrsBp1,
fisBS5-[:REGULATE {evidence:"Binding of purified proteins",evType:"strong",db:"RegulonDB",function:"activator"}]->rrsBp1,
fisBS6-[:REGULATE {evidence:"Binding of purified proteins",evType:"strong",db:"RegulonDB",function:"activator"}]->rrsBp1,
hnsBS1-[:REGULATE {evidence:"Binding of purified proteins",evType:"strong",db:"RegulonDB",function:"repressor"}]->rrsBp1,
hnsBS2-[:REGULATE {evidence:"Binding of purified proteins",evType:"strong",db:"RegulonDB",function:"repressor"}]->rrsBp1,
hnsBS3-[:REGULATE {evidence:"Binding of purified proteins",evType:"strong",db:"RegulonDB",function:"repressor"}]->rrsBp1,
fisBS1-[:NEXT]->rrsCp1,
fisBS2-[:NEXT]->rrsCp1,
fisBS3-[:OVERLAP]->rrsCp1,
fisBS4-[:NEXT]->rrsBp1,
fisBS5-[:NEXT]->rrsBp1,
fisBS6-[:OVERLAP]->rrsBp1,
hnsBS1-[:NEXT]->rrsBp1,
hnsBS2-[:OVERLAP]->rrsBp1,
hnsBS3-[:OVERLAP]->rrsBp1,
Fis-[:BINDS {evidence:"Binding of purified proteins",evType:"strong",db:"RegulonDB"}]->fisBS1,
Fis-[:BINDS {evidence:"Binding of purified proteins",evType:"strong",db:"RegulonDB"}]->fisBS2,
Fis-[:BINDS {evidence:"Binding of purified proteins",evType:"strong",db:"RegulonDB"}]->fisBS3,
Fis-[:BINDS {evidence:"Binding of purified proteins",evType:"strong",db:"RegulonDB"}]->fisBS4,
Fis-[:BINDS {evidence:"Binding of purified proteins",evType:"strong",db:"RegulonDB"}]->fisBS5,
Fis-[:BINDS {evidence:"Binding of purified proteins",evType:"strong",db:"RegulonDB"}]->fisBS6,
hns-[:BINDS {evidence:"Binding of purified proteins",evType:"strong",db:"RegulonDB"}]->hnsBS1,
hns-[:BINDS {evidence:"Binding of purified proteins",evType:"strong",db:"RegulonDB"}]->hnsBS2,
hns-[:BINDS {evidence:"Binding of purified proteins",evType:"strong",db:"RegulonDB"}]->hnsBS3;

match (tax:Taxon:Strain)-[:ORGANISM]->(chrom)
where tax.name="Escherichia coli str. K-12 substr. MG1655"
with chrom
match (feat:TFBS)
create chrom-[:CONTAINS]->feat;

