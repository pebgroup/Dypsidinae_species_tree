#!/usr/bin/env python

import phylo3, newick3, sys, os

if __name__ == "__main__":

    if len(sys.argv) < 5:
        print "usage: calc_branchwise_decisiveness.py <labeled_tree> <autophy_metadata_file> <node_scores_outfile> <locus_scores_outfile> <branch_counts_outfile>"
        sys.exit(0)

    # get infile names
    tree_file_path = sys.argv[1]
    autophy_metadata_file_path = sys.argv[2]

    # load tree
    with open(tree_file_path,"r") as treefile:
        tree = newick3.parse(treefile.readline())

    # extract locus sampling info from metadata file
    loci = {}
    with open(autophy_metadata_file_path,"r") as samplingfile:
        curlocus = ""
        first = True
        for line in samplingfile:
            if first:
                first = False
                print "assuming first line of locus file is column labels, it will be skipped"
                continue

            if len(line.strip()) == 0:
                continue

            line = line.strip()
            parts = line.split(",")
            tipname = parts[0]
            locusname = parts[1]
            if locusname not in loci.keys():
                loci[locusname] = []
            
            loci[locusname].append(tipname)
        
    nbrs = 0

    print """finding the loci that could inform the parent branch of each node"""
    with open(sys.argv[3],"w") as nodes_outfile:
        nodes_outfile.write("node_label,number_dec_loci"+",".join(loci.keys())+"\n")
        
        dec_loci_for_branches = {}
        n_brs_decisive_for_locus = {}
        for node in tree.iternodes():

            print node.label
            node.dec_loci = set()

            if len(node.children) > 2:
                print "node " + str(node) + " " + str(node.label) + " has more than 2 children!"
                continue

            if not node.istip and node.parent != None:
                left_child = node.children[0]
                right_child = node.children[1]

                for n in node.parent.children:
                    if n != node:
                        sister = n

                assert(sister)

            elif node.parent == None: # root of tree
                continue

            # this is brute force method. there should be a more efficient way
            # using recursion from tips to keep track of branches whose child branches
            # we've already validated

            for locus_name, locus_exemplars in loci.iteritems():

                if node.istip:
                    if node.label in locus_exemplars:
                        node.dec_loci.add(locus_name)
                    continue

                # have to have one exemplar of each child and one child of sister.
                # upstream exemplar can be provided by root (always known, so we don't check for it)
                left_ok = False
                right_ok = False
                sister_ok = False

                # check left
                for lchild in left_child.leaves():        
                    if lchild.label in locus_exemplars:
                        left_ok = True
                        break

                if not left_ok: # not decisive
                    continue

                else: # check right
                    for rchild in right_child.leaves():
                        if rchild.label in locus_exemplars:
                            right_ok = True
                            break

                if not right_ok: # not decisive
                    continue

                else: # check sister
                    for schild in sister.leaves():
                        if schild.label in locus_exemplars:
                            sister_ok = True
                            break

                if sister_ok: # decisive!
                    node.dec_loci.add(locus_name)
                    
                    if locus_name not in n_brs_decisive_for_locus:
                        n_brs_decisive_for_locus[locus_name] = 0
                    
                    n_brs_decisive_for_locus[locus_name] += 1
            
            # record info for this node
            nodes_outfile.write(",".join([node.label,str(len(node.dec_loci))]+[("1" if l in node.dec_loci else "0") for l in loci.keys()])+"\n")
            nbrs += 1

    # generate locus counts outfile
    with open(sys.argv[4],"w") as counts_outfile:
        counts_outfile.write("locus,prop_brs_decisive,num_brs_decisive\n")

        outdir = "individual_locus_decisiveness_trees/"
        try:
            os.mkdir(outdir)
        except OSError:
            pass

        for locusname, nbrs_decisive in n_brs_decisive_for_locus.iteritems():
            counts_outfile.write(",".join([locusname, str(float(nbrs_decisive)/nbrs), str(nbrs_decisive)]) + "\n")

    # generate branch counts output
    branch_counts = {}
    for n in tree.iternodes():
        ndecloci = len(n.dec_loci)
        if ndecloci == 0:
            continue
        
        if ndecloci not in branch_counts.keys():
            branch_counts[ndecloci] = 0

        branch_counts[ndecloci] += 1

    with open(sys.argv[5],"w") as branch_counts_outfile:
        branch_counts_outfile.write("nloci_decisive,n_branches\n")
        for n, c in branch_counts.iteritems():
            branch_counts_outfile.write(str(n) + "," + str(c) + "\n")
        branch_counts_outfile.close()
