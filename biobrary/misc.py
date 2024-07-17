import copy
import sys

def merge_islands(islands):
    """
    Merge islands

    Parameters
    --------------
    islands : a list of number contain islands need to merge
        A list in form like [[1, 3], [2, 4], [5, 9]]

    Returns
    ----------
    islands_merged : a list of number contain merged islands
    """
    
    labeled_border = []

    index = 0
    for isol in islands:
        labeled_border.append((isol[0], "l", index))
        labeled_border.append((isol[1], "r", index))
        index += 1
    labeled_border.sort(key=lambda x:x[0])
    left = 0
    right = 0
    islands_merged = []
    one_isol = []
    for border in labeled_border:
        if border[1] == "l":
            left += 1
        elif border[1] == "r":
            right += 1
        if left - right == 1 and len(one_isol) == 0:
            one_isol.append(border[0])
        elif left == right:
            one_isol.append(border[0])
            islands_merged.append(one_isol)
            one_isol = []

    return islands_merged


def island_location_by_position(islands: "list", pos: "int"):
    """
    Given a list of islands, and a point with absolut postion. Calculate
    which lsland this point located and the distances from both islands
    at both side.
    The position idx is start by 1.

    Paramters
    --------------
    islands: a sorted list of islands.
    pos: a int represent the absolute position

    Returns
    --------------
    location: the id of island, on which the point located.
        The id of left most island is 1. None for out of islands. 
    dist_left: the distance or indix from the island at most left side.
        None if location == None
    dist_right: the distance of index from the island at most right side.
        None if location  == None
    """
    location = None
    idx = 0
    for island in islands:
        idx += 1
        if pos >= island[0] and pos <= island[1]:
            location = idx

    if not location:
        idx_left = None
        idx_right = None
    else:
        island_size = []
        for island in islands:
            island_size.append(island[1] - island[0] + 1)
        idx_left = sum(island_size[:location - 1]) + (pos - islands[location - 1][0] + 1)
        idx_right = sum(island_size[location:]) + (islands[location - 1][1] - pos + 1)
    
    return location, idx_left, idx_right


def island_location_by_distance(islands: "list", dist: "int", side="left"):
    """
    Given a list of lslands, and a distance of a point form left side of
    left most island or a distance from right side of right most island.
    Calculate the location of this point, and its relative postion to the
    islands

    Parameters
    --------------
    islands: a sorted list of islands.
    dist: a int represent the offset from the side.
    side: which lsland the distance relatived to.
    
    Returns
    ---------------
    location: the id of island, on which the point located.
        The id of left most island is 1. None for out of islands. 
    pos: relative position rative to the side.
    """
    location = None
    if side == "left":
        cumsum = [islands[0][1] - islands[0][0] + 1]
        for idx, island in enumerate(islands[1:]):
            cumsum.append(cumsum[idx] + island[1] - island[0] + 1)

        for idx, cum in enumerate(cumsum):
            if cum > dist:
                location = idx + 1
                break
            if location:
                pos = islands[location - 1][0] + (pos - cumsum[location - 2])

    elif side == "right":
        cumsum = [islands[-1][1] - islands[-1][0] + 1]
        for idx, island in enumerate(islands[:-1][::-1]):
            cumsum.append(cumsum[idx] + island[1] - island[0] + 1)

        for idx, cum in enumerate(cumsum):
            if cum > dist:
                location = idx + 1
                break
            if location:
                pos = islands[-location][1] - (pos - cumsum[location - 2]) 
    else:
        print(f"{side} not is not a acceptable parameter.", file=sys.stderr)
    
    if not location:
        pos = None
    
    return location, pos


def island_distance(islands: "list", location: "int", side="left") -> "int":
    """
    Calculate distance for a given posistion.

    Parameters
    --------------
    islands: list.
        A sorted list of integer. 'left <= right' for every lands.
    location: int
    side: left or right
        The distance from the left side or right side.

    Retures
    ------------
    distance: int or None
        None for the location not on lsolands

    """
    dist = 0
    if side == "left":
        for isol in islands:
            if location >= isol[0] and location <= isol[1]:
                dist += location - isol[0] + 1
                return dist
            elif location > isol[1]:
                dist += isol[1] - isol[0] + 1
            else:
                return None

    elif side == "right":
        for isol in islands[::-1]:
            if location >= isol[0] and location <= isol[1]:
                dist += isol[1] - location + 1
                return dist
            elif location < isol[0]:
                dist += dist[1] - isol[0] + 1
            else:
                return None
    else:
        print(f"{side} not know.")
        return 0


def island_location(islands: "list", distance: "int", side="left") -> "int":
    """
    Get the location of a point by distance.

    Parameters
    ------------
    islands: list
        A sorted list of int, 'left <= right' for every island.
    distance: int
        Distance from the left or right.

    Retures
    -------------
    location: int
        Location on the islands.
    """
    if side == "left":
        for isol in islands:
            distance -= isol[1] - isol[0] + 1
            if distance <= 0:
                return distance + isol[1]
    elif side == "right":
        for isol in islands[::-1]:
            distance -= isol[1] - isol[0] + 1
            if distance <= 0:
                return isol[0] - distance
    else:
        print(f"{side} not know.")
        return 0


def nucleotide_relative_coord(positions: "list", reference: "int", reverse: "bool" =False) -> "list":
    """
    Calculate relavtive location to the reference.

    Parameters
    -----------
    reference: int.
        To which all other locatation ralative to.
    positions: list.
    reverse: bool, default is False.
        If set to True, the sign of location would reversed.

    Retures
    -----------
    relative_pos: list.
    """

    if reverse:
        pos_rela = [reference - pos for pos in positions]
    else:
        pos_rela = [pos - reference for pos in positions]
    pos_rela = [ele + 1 if ele > -1 else ele for ele in pos_rela]
    return pos_rela


def change_coordinate_rel(pos, left, right, ori="+"):
    """
    Convert base cooradition

    Parameters:
    ------------------
    pos: list, the pos to be converted.[[l1, l2], [l3, l4],...]
    left: left side
    right: right side
    ori: "+"(default) or "-".
        oritation

    Returns
    --------------
    list: converted result.
    """
    if ori == '+':
        pos.sort(key=lambda x:x[0])
        pos_n = []
        for p in pos:
            pos_n += p
        pos_n = nucleotide_relative_coord(pos_n, left, reverse=False)
        pos = []
        for idx in range(0, len(pos_n), 2):
            pos.append([pos_n[idx], pos_n[idx + 1]])
        return pos
    else:
        pos.sort(key=lambda x:x[1])
        pos_n = []
        for p in pos:
            pos_n += p
        pos_n = nucleotide_relative_coord(pos_n, right, reverse=True)
        pos_n.reverse()
        pos = []
        for idx in range(0, len(pos_n), 2):
            pos.append([pos_n[idx], pos_n[idx + 1]])
        return pos


def change_coordinate_abs(ref, pos, ori="+"):
    ref = copy.deepcopy(ref)
    ref.sort(key=lambda x:x[0])
    if ori == "+":
        for ele in ref:
            block_len = ele[1] - ele[0] + 1
            if pos <= block_len:
                return ele[0] + pos - 1
            pos -= block_len

    else:
        for ele in ref[::-1]:
            block_len = ele[1] - ele[0] + 1
            if pos <= block_len:
                return ele[1] - pos + 1
            pos -= block_len


def exon_coordinate_rel(exon, pos):
    mark = 0
    pos_exon = 0
    for exon_ele in exon:
        if pos >= exon_ele[0] and pos <= exon_ele[1]:
            pos_exon += pos - exon_ele[0] + 1
            mark = 1
            break
        else:
            pos_exon += exon_ele[1] - exon_ele[0] + 1
    if mark:
        return pos_exon
    else:
        return None


def split_block(segments, block_start, block_len):
    data_out = []
    block_len_left = block_len
    seg_len = len(segments)
    start_in = 0
    for idx in range(seg_len):
        if segments[idx][0] <= block_start and segments[idx][1] >= block_start:
            start_in = 1
            break
    if not start_in:
        print("Waring, start not in segments, [] returned.")
        return []

    while True:
        consum_len = segments[idx][1] - block_start + 1
        block_len_left -= consum_len
        if block_len_left > 0:
            if idx < seg_len - 1:
                data_out.append([block_start, segments[idx][1]])
                idx += 1
                block_start = segments[idx][0]
            else:
                data_out.append([block_start, segments[idx][1]])
                break
        else:
            data_out.append([block_start, block_start + block_len - 1])
            break
        block_len = block_len_left
    if block_len_left > 0:
        print("Warning, block out of segments.")

    return data_out


if __name__ == "__main__":
    print(change_coordinate_rel([[21, 30], [42, 60], [79, 85]], 21, 103, ori='+'))
    print(change_coordinate_rel([[21, 30], [42, 60], [79, 85]], 21, 103, ori='-'))
