input_file = open('/home/artem/work/2014/PhysRev/doc_list.txt', 'r')
#output_file = open('out.bib', 'w')
lines = input_file.readlines()
i = 1
err = 0
for line in lines:
    line = line.split('] ')[1]
    if 'and' in line:
        spl_authors = line.split(' and ')
        spl_other = spl_authors[1].split(',')
        authors = spl_authors[0] + ' and ' + spl_other[0]
        if ', and' in authors:
            authors = authors.replace(', and', ' and')
        alias = authors.split(',')[0].split(' ')[-1]
    else:
        authors = line.split(',')[0]
        alias = authors.split(' ')[-1]
    #print alias, '\n', authors
    print spl_other
    print i
    i += 1
input_file.close()
