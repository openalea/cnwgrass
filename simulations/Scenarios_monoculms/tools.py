import os
import pandas as pd
from collections import deque


def is_int(input):
    try:
        num = int(input)
    except ValueError:
        return False
    return True


def _buildDic(keyList, val, dic):
    if len(keyList) == 1:
        dic[keyList[0]] = val
        return

    newDic = dic.get(keyList[0], {})
    dic[keyList[0]] = newDic
    _buildDic(keyList[1:], val, newDic)


def buildDic(dict_scenario, dic=None):
    if not dic:
        dic = {}

    for k, v in dict_scenario.items():
        if not pd.isnull(v):
            keyList = k.split(':')
            keyList_converted = []
            for kk in keyList:
                if is_int(kk):
                    keyList_converted.append(int(kk))
                else:
                    keyList_converted.append(kk)
            _buildDic(keyList_converted, v, dic)

    return dic


def nested_enumerate(lst):
    """An analogue of enumerate for nested lists.

       Returns an iterator over the (index, element) pairs of `lst` where
       `index` is a list of integers [x0, x1,.., xn] such that
       `lst[x0][x1]...[xn]==element`


       >>> for i, e in nested_enumerate([0, [[[1, [2, [[[3]]]]]]], [[[4]]]]):
               print('%s %s'%(str(i), str(e)))
       [0] 0
       [1, 0, 0, 0] 1
       [1, 0, 0, 1, 0] 2
       [1, 0, 0, 1, 1, 0, 0, 0] 3
       [2, 0, 0, 0] 4
    """

    # initial, partial index of lst
    partial_index = deque([([i], e) for (i, e) in enumerate(lst)])

    while partial_index:
        index, obj = partial_index.popleft()
        if isinstance(obj, list):
            # if obj is a list then its elements require further indexing
            new_dimension = [(index + [i], e) for (i, e) in enumerate(obj)]
            partial_index.extendleft(reversed(new_dimension))
        else:
            # obj is fully indexed
            yield index, obj


def clear_directory(directory_path):
    for file_name in os.listdir(directory_path):
        file_path = os.path.join(directory_path, file_name)
    try:
        if os.path.isfile(file_path):
            os.unlink(file_path)
    except Exception as e:
        print(e)
