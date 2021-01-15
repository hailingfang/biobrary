
import ete3


class CircleNode:
    import sys

    def __init__(self, BaseInnerNode, Name=''):
        self.name = Name
        self.attr = None
        self.base_dist = None
        self.linker_dist = None
        #self.base_inner_node can be 'None' or a ete3 TreeNode
        self.base_inner_node = BaseInnerNode
        #self.birth_inner_node can be '[]' or a list of ete3 TreeNode
        self.birth_inner_node = set()
        # self.base_link_bridge can be 'None' or '[BaseInnerNode, parent_circle_node]'
        self.base_link_bridge = None
        # self.birth_link_bridge can be '[]' or '[[birth_inner_node, child_circle_node], ...]'
        self.birth_link_bridge = set()

    def is_leaf(self):
        if self.birth_inner_node:
            return False
        else:
            return True

    def is_child(self, second_circle_node):
        children = self.get_children()
        if second_circle_node in children:
            return True
        else:
            return False

    def is_parent(self, second_circle_node):
        parent = self.get_parent()
        if second_circle_node == parent:
            return True
        else:
            return False

    def is_root(self):
        if self.base_link_bridge:
            return False
        else:
            return True

    def is_link_directly(self, second_circle_node):
        all_linked_circle_node = self.birth_inner_node.append(self.parent_link_bridge)
        if second_circle_node in all_linked_circle_node:
            return True
        else:
            return False

    def is_base_inner_node(self, inner_node):
        if inner_node == self.base_inner_node:
            return True
        else:
            return False

    def is_birth_inner_node(self, inner_node):
        if inner_node in self.birth_inner_node:
            return True
        else:
            return False

    def get_children(self):
        children = []
        for bridge in self.birth_link_bridge:
            children.append(bridge[1])
        return children

    def get_parent(self):
        if self.base_link_bridge:
            return self.base_link_bridge[1]
        else:
            return None

    def get_children_base_inner_node(self):
        children = self.get_children()
        children_base_inner_node = [child.base_inner_node for child in children]
        return children_base_inner_node

    def get_inner_node(self):

        def get_inner_node_recursive(parent, stop_node):
            dt_out = []
            if parent not in stop_node:
                dt_out.append(parent)
                for child in parent.children:
                    dt_out += get_inner_node_recursive(child, stop_node)
            return dt_out

        base_inner_node = self.base_inner_node
        children_base_inner_node = self.get_children_base_inner_node()
        inner_node = get_inner_node_recursive(base_inner_node, children_base_inner_node)
        return inner_node

    def get_root_circle_node(self):
        root = self
        parent_link = root.base_link_bridge
        while parent_link:
            root = parent_link[1]
            parent_link = root.base_link_bridge
        return root

    def get_birth_node_to_child(self, second_circle_node):
        children = self.get_children()
        if second_circle_node in children:
            birth_node = [bridge[0] for bridge in self.birth_link_bridge if \
                bridge[1] == second_circle_node][0]
            return birth_node
        else:
            return None

    def get_parent_birth_node(self):
        parent = self.get_parent()
        birth_node = parent.get_birth_node_to_child(self)
        return birth_node

    def get_descendants(self):
        descendants = []
        for child in self.get_children():
            descendants.append(child)
            descendants += child.get_descendants()
        return descendants

    def add_child(self, parent_inner_link_node, child_circle_node):
        self.birth_link_bridge.add((parent_inner_link_node, child_circle_node))
        child_circle_node.base_link_bridge = (child_circle_node.base_inner_node, self)

    def link_to_parent(self, parent_inner_link_node, parent_circle_node):
        self.base_link_bridge = [self.base_inner_node, parent_circle_node]
        parent_circle_node.birth_link_bridge.append((parent_inner_link_node, self))

    def detach(self):
        parent = self.get_parent()
        if parent:
            for bridge in parent.birth_link_bridge:
                if bridge[1] == self:
                    tmp = bridge
                    break
            parent.birth_link_bridge.remove(tmp)
        return self

    def delete(self):
        children = self.get_children()
        parent = self.get_parent()
        parent_birth_node = self.get_parent_birth_node()
        if not self.is_root():
            for child in children:
                child.base_dist += self.base_dist
                child.linker_dist = self.linker_dist + child.base_dist
                parent.add_child(parent_birth_node, child)
            self.detach()
        return self

    def merge(self, second_circle_node):
        if self.is_child(second_circle_node):
            birth_node = self.get_birth_node_to_child(second_circle_node)
            self.birth_inner_node.remove(birth_node)
            self.birth_link_bridge += second_circle_node.birth_link_bridge
            second_circle_node.detach()
        elif self.is_parent(second_circle_node):
            birth_node = self.get_parent_birth_node(second_circle_node)
            second_circle_node.birth_inner_node.remove(birth_node)
            second_circle_node.birth_link_bridge += self.birth_link_bridge
            self.detach()

    def make_profile_tree(self):

        def traverse_tree(circle_node, parent_node):
            for child in circle_node.get_children():
                child_node = ete3.TreeNode()
                child_node.name = child.name
                child_node.dist = child.base_dist
                parent_node.add_child(child_node)
                traverse_tree(child, child_node)
            return 0

        root_circle_node = self.get_root_circle_node()
        root = ete3.TreeNode()
        root.dist = 0
        root.name = root_circle_node.name
        traverse_tree(root_circle_node, root)
        return root

    def split_inner_tree(self):
        root_circle_node = self.get_root_circle_node()
        descendants = root_circle_node.get_descendants()
        for des in descendants:
            base_inner_node = des.base_inner_node
            base_inner_node.detach()
        return 0

    def join_inner_tree(self):
        root_circle_node = self.get_root_circle_node()
        descendants = root_circle_node.get_descendants()
        for node in descendants:
            base_inner_node = node.base_inner_node
            birth_node = node.get_parent_birth_node()
            birth_node.add_child(base_inner_node)
        return 0

    def trim_inner_node(self):

        def trim_node(node, stop_node):
            nonlocal child_need_detach
            nonlocal base_inner_node
            if node not in stop_node:
                i = 0
                for node_ele in node.get_descendants():
                    if node_ele in stop_node:
                        i = 1
                        break
                if (i == 0) and (node != base_inner_node):
                    child_need_detach.append(node)
                else:
                    for child in node.children:
                        trim_node(child, stop_node)
            return child_need_detach

        base_inner_node = self.base_inner_node
        stop_node = self.get_children_base_inner_node()
        child_need_detach = []
        trim_node(base_inner_node, stop_node)
        for node in child_need_detach:
            node.detach()

    def traverse(self):
        nodes = []
        nodes.append(self)
        for child in self.get_children():
            nodes += child.traverse()
        return nodes

    def change_base_inner_node_name(self):
        root_circle_node = self.get_root_circle_node()
        for node in root_circle_node.traverse():
            if not node.base_inner_node.is_leaf():
                node.base_inner_node.name = node.name
        return 0

    def print_cluster(self, f_out=sys.stdout):
        if f_out != sys.stdout:
            f_out = open(f_out, 'w')
        root_circle_node = self.get_root_circle_node()
        root_circle_node.split_inner_tree()
        for node in root_circle_node.traverse():
            print('>' + node.name, file=f_out)
            print('@' + node.base_inner_node.write(), file=f_out)
            inner_node = node.get_inner_node()
            for inner_node_ele in inner_node:
                if inner_node_ele.name != '':
                    print(inner_node_ele.name, file=f_out)
        root_circle_node.join_inner_tree()
        f_out.close()
        return 0


class CircleNodeTree:

    def __init__(self, treefile, edge_len_cutoff):

        def find_birth_linker(parent_node, edge_len_cutoff):
            dt_out = []
            for child in parent_node.children:
                if child.dist > edge_len_cutoff:
                    dt_out.append([parent_node, child])
                else:
                    sub_dt = find_birth_linker(child, edge_len_cutoff)
                    dt_out += sub_dt
            return dt_out

        def link_circle_node(parent_circle_node, edge_len_cutoff):
            nonlocal i
            birth_linker = find_birth_linker(parent_circle_node.base_inner_node, edge_len_cutoff)
            for linker in birth_linker:
                i += 1
                parent_inner_link_node = linker[0]
                child_base_inner_node = linker[1]
                child_circle_node_name = 'CircleNode_' + str(i)
                child_circle_node = CircleNode(child_base_inner_node, child_circle_node_name)
                child_circle_node.name = child_circle_node_name
                child_circle_node.base_dist = child_base_inner_node.get_distance(parent_circle_node.base_inner_node)
                child_circle_node.linker_dist = child_base_inner_node.get_distance(parent_inner_link_node)
                parent_circle_node.birth_inner_node.add(parent_inner_link_node)
                parent_circle_node.add_child(parent_inner_link_node, child_circle_node)
                link_circle_node(child_circle_node, edge_len_cutoff)

        self.original_tree = ete3.Tree(treefile)
        i = 0
        root_circle_node = CircleNode(self.original_tree, 'CircleNode_0')
        link_circle_node(root_circle_node, edge_len_cutoff)
        self.circle_node_tree = root_circle_node
