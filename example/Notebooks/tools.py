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
