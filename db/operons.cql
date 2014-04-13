CREATE (gltT:Gene:Feature {name:"gltT",blattner:"b3969"}),(rrfB:Gene:Feature {name:"rrfB",blattner:"b3971"}),
(rrlB:Gene:Feature {name:"rrlB",blattner:"b3970"}),(rrsB:Gene:Feature {name:"rrsB",blattner:"b3968"}),
(gltU:Gene:Feature {name:"gltU",blattner:"b3757"}),(rrfC:Gene:Feature {name:"rrfC",blattner:"b3759"}),
(rrlC:Gene:Feature {name:"rrlC",blattner:"b3758"}),(rrsC:Gene:Feature {name:"rrsC",blattner:"b3756"}),
(rrsBp1:Promoter:Feature {name:"rrsBp1"}),(rrsBp2:Promoter:Feature {name:"rrsBp2"}),
(rrsCp1:Promoter:Feature {name:"rrsCp1"}),(rrsCp2:Promoter:Feature {name:"rrsCp2"}),
(op1:Operon {name:"rrsB-gltT-rrlB-rrfB"}),(op2:Operon {name:"rrsC-gltU-rrlC-rrfC"}),
(tu1:TUnit {name:"rrsB-gltT-rrlB-rrfB"}),(tu2:TUnit {name:"rrsC-gltU-rrlC-rrfC"}),
rrsBp1-[:CONTROL {evidence:"Boundaries of transcription experimentally identified",evType:"weak",db:"RegulonDB"}]->tu1,rrsBp2-[:CONTROL]->tu1,
rrsCp1-[:CONTROL {evidence:"Boundaries of transcription experimentally identified",evType:"weak",db:"RegulonDB"}]->tu2,rrsCp2-[:CONTROL]->tu2,
gltT-[:PART_OF]->tu1,rrfB-[:PART_OF]->tu1,rrlB-[:PART_OF]->tu1,rrsB-[:PART_OF]->tu1,
gltU-[:PART_OF]->tu2,rrfC-[:PART_OF]->tu2,rrlC-[:PART_OF]->tu2,rrsC-[:PART_OF]->tu2,
tu1-[:PART_OF]->op1,tu2-[:PART_OF]->op2
;

match (tax:Taxon:Strain)-[:ORGANISM]->(chrom)
where tax.name="Escherichia coli str. K-12 substr. MG1655"
with chrom
match (feat:Feature)
create chrom-[:CONTAINS]->feat;

