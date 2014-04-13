CREATE (troot:Taxon { id:131567, name:"cellular organisms", rank:"root"}),
(sking:Taxon { id:2, name:"Bacteria", rank:"superkingdom"}),
(phyl:Taxon { id:1224, name:"Proteobacteria", rank:"phylum"}),
(cls:Taxon { id:1236, name:"Gammaproteobacteria", rank:"class"}),
(ord:Taxon { id:91347, name:"Enterobacteriales", rank:"order"}),
(fam:Taxon { id:543, name:"Enterobacteriaceae", rank:"family"}),
(gen:Taxon:Genus { id:561, name:"Escherichia", rank:"genus"}),
(spc:Taxon:Species { id:562, name:"Escherichia coli", rank:"species"}),
(subsp:Taxon { id:83333, name:"Escherichia coli K-12"}),
(str:Taxon:Strain { id:511145, name:"Escherichia coli str. K-12 substr. MG1655"}),
str-[:IS_A]->subsp,subsp-[:IS_A]->spc,spc-[:IS_A]->gen,gen-[:IS_A]->fam,fam-[:IS_A]->ord,ord-[:IS_A]->cls,cls-[:IS_A]->phyl,phyl-[:IS_A]->sking,sking-[:IS_A]->troot;
CREATE INDEX ON :Taxon(name);

CREATE (chrom:Genome {id:"NC_000913.2", name:"Escherichia coli str. K-12 substr. MG1655, complete genome", length:4639675, geom:"circular", timestmp:"22-JUL-2013"})
with chrom
match (tax:Taxon:Strain)
where tax.name="Escherichia coli str. K-12 substr. MG1655"
CREATE (tax)-[:ORGANISM]->(chrom);
