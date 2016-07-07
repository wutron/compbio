# see also compbio.phylo

# TODO
# - transfers only make sense for binary trees
#   ==> generalize to non-binary trees
# - currently the recipient species for a transfer edge is ALWAYS equal to
#   the reconciled species for the child of the transfer edge
#   (this is okay for non-dated species trees because this is most parsimonious)
#   ==> generalize to separate transfer edges from donor and recipient species
#   add trans_sp (similar to trans) that contains the recipient species for a transfer edge
#   let trans_sp default to trans for above case (recipient species for transfer edge ==
#     reconciled species for child of transfer edge)

from rasmus import util

#=============================================================================
# Reconciliation functions

def find_loss_node(node, recon, events, trans):
    """Finds the loss events for a branch in a reconciled gene tree"""
    loss = []

    # if not parent, then no losses
    if not node.parent:
        return loss

    # determine starting and ending species
    sstart = recon[node]
    # special case if transfer edge
    # TODO: currently never infers any losses since sstart == send
    #       need to look at recon[trans_sp[node.parent]]i
    if (node.parent in trans) and (trans[node.parent] == node):
        send = recon[trans[node.parent]]
    else:
        send = recon[node.parent]

    # determine species path of this gene branch (node, node.parent)
    ptr = sstart
    spath = []
    while ptr != send:
        spath.append(ptr)
        ptr = ptr.parent

    # determine whether node.parent is a dup or a trans
    # if so, send (species end) is part of species path
    if events[node.parent] == "dup" or events[node.parent] == "trans":
        spath.append(send)

    # go up species path (skip starting species)
    # every node on the list is at least one loss
    for i, snode in enumerate(spath[1:]):
        for schild in snode.children:
            if schild != spath[i]:
                loss.append([node, schild])

    return loss

def find_loss(gtree, recon, events, trans):
    """Returns a list of gene losses in a gene tree"""
    loss = []

    def walk(node):
        loss.extend(find_loss_node(node, recon, events, trans))
        node.recurse(walk)
    walk(gtree.root)

    return loss

def count_dup(gtree, events):
    """Returns the number of duplications in a gene tree"""
    ndups = 0

    def walk(node):
        if events[node] == "dup":
            ndups += len(node.children) - 1
        node.recurse(walk)
    walk(gtree.root)

    return ndups

def count_trans(gtree, events):
    """Returns the number of transfers in a gene tree"""
    ntrans = 0

    def walk(node):
        if events[trans] == "trans":
            assert len(node.children) == 2
            ntrans += 1
        node.recurse(walk)
    walk(gtree.root)

    return ntrans

def count_loss(gtree, recon, events, trans):
    """Returns the number of losses in a gene tree"""
    return len(find_loss(gtree, recon, events, trans))

def count_dup_trans_loss(gtree, recon, events, trans):
    """Returns the number of duplications + transfers + losses in a gene tree"""
    ndups = count_dup(gtree, events)
    ntrans = count_trans(gtree, events)
    nloss = count_loss(gtre, recon, events, trans)
    return ndups + ntrans + nloss

def find_xenologs(gtree, stree, recon, events, trans, counts=True, species_branch=False):
    """Find all xenolog pairs within a gene tree

    NOTE: THIS HAS NOT BEEN TESTED!!!
    """
    xenos = []

    for node, event in events.items():
        if event == "trans":
            assert len(node.children) == 2
            if trans[node] == node.children[0]:
                children = (node.children[1], node.children[0])
            else:
                children = node.children
            leavesmat = [x.leaves() for x in children]
            sp_counts = [util.hist_dict(util.mget(recon, row))
                         for row in leavesmat]

            for i in range(len(leavesmat)):
                for j in range(i+1, len(leavesmat)):
                    for gene1 in leavesmat[i]:
                        for gene2 in leavesmat[j]:
                            g1, g2 = gene1, gene2
                            a, b = i, j

                            xeno = [g1.name, g2.name]
                            if counts:
                                xeno.extend([sp_counts[a][recon[g1]],
                                             sp_counts[b][recon[g2]]])
                            if species_branch:
                                xeno.append(recon[node])
                            xenos.append(tuple(xenos))

    return xenos

def subset_recon(tree, recon, events=None, trans=None):
    """Ensure the reconciliation only refers to nodes in tree"""

    # get all nodes that are walkable
    nodes = set(tree.postorder())

    for node in list(recon):
        if node not in nodes:
            del recon[node]

    if events:
        for node in list(events):
            if node not in nodes:
                del events[node]

    if trans:
        for node in list(trans):
            if node not in nodes:
                del trans[node]

        for node, tnode in trans.iteritems():
            # if transfer edge has been lost,
            # walk down until we find a node that is still in the tree
            while tnode not in nodes:
                assert len(tnode.children) == 1
                tnode = tnode.children[0]
            trans[node] = tnode


#=============================================================================
# Reconciliation Input/Output

def assert_events(events, trans):
    assert all([node.is_leaf() for (node, event) in events.iteritems() if event == "gene"])
    assert set(trans) == set([node for (node, event) in events.iteritems() if event == "trans"])
    assert len(trans) == 0 or \
           all([trans[node].parent is node
                for (node, event) in events.iteritems() if event == "trans"])

def write_recon_events(filename, recon, events, trans):
    """Write a reconciliation and events to a file"""
    assert_events(events, trans)

    util.write_delim(filename,
                     [(str(node.name), str(snode.name), events[node])
                      + ((trans[node].name,) if node in trans else ())
                      for node,snode in recon.items()])

def read_recon_events(filename, tree, stree):
    """Read a reconciliation and events data structure from file"""
    recon = {}
    events = {}
    trans = {}  # child of transfer edge

    for toks in util.read_delim(filename):
        if len(toks) == 3:
            name, sname, event = toks
            tname = None
        elif len(toks) == 4:
            name, sname, event, tname = toks
        else:
            raise Exception("invalid line: %s" % '\t'.join(toks))

        if name.isdigit(): name = int(name)
        if sname.isdigit(): sname = int(sname)
        node = tree.nodes[name]
        try:
            recon[node] = stree[sname]
        except:
            recon[node] = stree[str(sname)]
        events[node] = event
        if tname:
            if tname.isdigit(): tname = int(tname)
            trans[node] = tree.nodes[tname]

    assert_events(events, trans)

    return recon, events, trans


#============================================================================
# duplication transfer loss counting
#

def init_dup_trans_loss_tree(stree):
    """initalize counts to zero"""

    def walk(node):
        node.data['dup'] = 0
        node.data['trans'] = dict.fromkeys(stree.nodes, 0)
        node.data['loss'] = 0
        node.data['appear'] = 0
        node.data['genes'] = 0
        node.recurse(walk)
    walk(stree.root)

def count_dup_trans_loss_tree(tree, stree, gene2species, recon, events, trans, losses=None):
    """count dup trans loss"""

    if losses is None:
        losses = find_loss(tree, stree, recon, events, trans)

    ndup = 0
    ntrans = 0
    nloss = 0
    nappear = 0

    # count appearance
    recon[tree.root].data["appear"] += 1
    nappear += 1

    # count dups trans genes
    for node, event in events.iteritems():
        if event == "dup":
            recon[node].data['dup'] += 1
            ndup += 1
        elif event == "trans":
            recon[node].data['trans'][recon[trans[node]].name] += 1
            ntrans += 1
        elif event == "gene":
            recon[node].data['genes'] += 1

    # count losses
    for gnode, snode in losses:
        snode.data['loss'] += 1
        nloss += 1

    return ndup, ntrans, nloss, nappear

def count_ancestral_genes(stree):
    """count ancestral genes"""
    # TODO: figure out how to incorporate transfers into the count
    raise Exception("not implemented")
