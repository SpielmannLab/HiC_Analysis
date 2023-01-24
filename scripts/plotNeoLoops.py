#!/bin/python

############################
### Code created by Xiaotao Wang and James Hawley and published at
### https://github.com/XiaoTaoWang/NeoLoopFinder
### Modified by K. Schultz, January 2023
############################


from neoloop.visualize.core import * 
import cooler
import sys
clr = cooler.Cooler(sys.argv[1])
ids=[]
List = [line.rstrip() for line in open('allOnco-genes.txt')]
with open('%s'%sys.argv[2]) as file:
	for line in file:
		line=(line.rstrip())
		idstring=line.split('\t')[6] 
		id=idstring.partition(',')[0]
		ids.append(id)
with open('%s'%sys.argv[3]) as file:
	for line in file:
		assembly=(line.rstrip())
		a_id=assembly.partition('\t')[0] 
		for id in ids :
			if id == a_id :
				vis = Triangle(clr, assembly, n_rows=5, figsize=(7, 5.2), track_partition=[5, 0.8, 0.8, 0.2, 0.5], correct='weight', span=300000, space=0.08)
				vis.matrix_plot(vmin=0, cbr_fontsize=9)
				vis.plot_chromosome_bounds(linewidth=2)
				#vis.plot_genes(release=75, filter_=List, fontsize=10)
				vis.plot_genes(release=75, fontsize=10)
				vis.plot_chromosome_bar(name_size=10, coord_size=5)
				vis.outfig('%s.%s.pdf'%(sys.argv[4], id))
				vis.plot_loops('%s'%sys.argv[2], face_color='none', marker_size=40, cluster=True)
				vis.outfig('%s.%s.marked.pdf'%(sys.argv[4], id))
				plt.close('all')
