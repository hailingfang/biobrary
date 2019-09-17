import ete3

def gen_test_tree():
    a = ete3.TreeNode()
    a.name = 'root'
    a.dist = 0

    b = ete3.TreeNode()
    b.name = 'b'
    b.dist = 2

    c = ete3.TreeNode()
    c.name = 'c'
    c.dist = 1

    a.add_child(b)
    a.add_child(c)

    d = ete3.TreeNode()
    d.name = 'd'
    d.dist = 2

    e = ete3.TreeNode()
    e.name = 'e'
    e.dist = 2

    b.add_child(d)
    b.add_child(e)

    f = ete3.TreeNode()
    f.name = 'f'
    f.dist = 1

    g = ete3.TreeNode()
    g.name = 'g'
    g.dist = 1

    c.add_child(f)
    c.add_child(g)

    h = ete3.TreeNode()
    h.name = 'h'
    h.dist = 1

    i = ete3.TreeNode()
    i.name = 'i'
    i.dist = 1

    d.add_child(h)
    d.add_child(i)

    j = ete3.TreeNode()
    j.name = 'j'
    j.dist = 1

    k = ete3.TreeNode()
    k.name = 'k'
    k.dist = 1

    h.add_child(j)
    h.add_child(k)

    l = ete3.TreeNode()
    l.name = 'l'
    l.dist = 1

    m = ete3.TreeNode()
    m.name = 'm'
    m.dist = 2

    e.add_child(l)
    e.add_child(m)

    n = ete3.TreeNode()
    n.name = 'n'
    n.dist = 1

    o =ete3.TreeNode()
    o.name = 'o'
    o.dist = 1

    m.add_child(n)
    m.add_child(o)

    p = ete3.TreeNode()
    p.name = 'p'
    p.dist = 1

    q = ete3.TreeNode()
    q.name = 'q'
    q.dist = 2

    f.add_child(p)
    f.add_child(q)

    r = ete3.TreeNode()
    r.name = 'r'
    r.dist = 1

    s = ete3.TreeNode()
    s.name = 's'
    s.dist = 1

    q.add_child(r)
    q.add_child(s)

    t = ete3.TreeNode()
    t.name = 't'
    t.dist = 1

    u = ete3.TreeNode()
    u.name = 'u'
    u.dist = 1

    s.add_child(t)
    s.add_child(u)

    #print(a.write(features=['name']))
    #print(a)
    #print(a.name)
    return a

test_tree = gen_test_tree()
