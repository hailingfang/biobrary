import copy

def merge_isolands(isolands):
    """
    Merge isolands

    Parameters
    --------------
    isolands : a list of number contain isolands need to merge
        A list in form like [[1, 3], [2, 4], [5, 9]]

    Returns
    ----------
    isolands_merged : a list of number contain merged isolands
    """
    
    labeled_border = []

    index = 0
    for isol in isolands:
        labeled_border.append((isol[0], "l", index))
        labeled_border.append((isol[1], "r", index))
        index += 1

    labeled_border += 1

    left = 0
    right = 0
    isolands_merged = []
    isolands_merged.append(labeled_border[0][0])
    one_isol = []
    for border in labeled_border:
        if border[1] == "l":
            left += 1
        elif border[1] == "r":
            right += 1
        if left - right == 1:
            one_isol.append(border[0])
        if left == right:
            one_isol.append(border[0])
            isolands_merged.append(one_isol)
            one_isol = []

    return isolands_merged


def change_coordinate_rel(ref, pos, ori="+"):
    """
    Convert base cooradition

    Parameters:
    ------------------
    ref: [left, right], reference position.
    pos: list, the pos to be converted.[[l1, l2], [l3, l4],...]
    ori: "+"(default) or "-".
        oritation

    Returns
    --------------
    list: converted result.
    """
    ref = copy.deepcopy(ref)
    if type(pos) == list and type(pos[0]) == list:
        new_pos = []
        if ori == "+":
            for ele in pos:
                ele.sort()
            pos.sort(key=lambda x:x[0])
        else:
            for ele in pos:
                ele.sort(reverse=True)
            pos.sort(key=lambda x:x[0], reverse=True)

        if ori == "+":
            for pos_ele in pos:
                tmp = []
                for num in pos_ele:
                    tmp.append(num - ref[0] + 1)
                new_pos.append(tmp)
        else:
            for pos_ele in pos:
                tmp = []
                for num in pos_ele:
                    tmp.append(ref[1] - num + 1)
                new_pos.append(tmp)
        return new_pos
    
    elif type(pos) == list and type(pos[0]) == int:
        new_pos = []
        if ori == "+":
            for num in pos:
                new_pos.append(num - ref[0] + 1)
        else:
            for num in pos:
                new_pos.append(ref[1] - num + 1)
        return new_pos

    elif type(pos) == int:
        new_pos = None
        if ori == "+":
            new_pos = pos - ref[0] + 1
        else:
            new_pos = ref[1] - pos + 1

        return new_pos

    else:
        print("pos type not correct.")


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
        print("start not in segments")
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
        print("block out of segments")

    return data_out


if __name__ == "__main__":
    print(change_coordinate([300,1000], [3,4,9], ori="-", coor="abs"))
    print(split_block([[5, 8],[10, 15],[20, 28]], 7, 3))
