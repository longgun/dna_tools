import glob
import timeit

import pytest

from dnase import DnaseTools


@pytest.fixture(scope="function")
def dnase_obj():
    res = DnaseTools()
    yield res
    del res


def test_time(dnase_obj):
    """time test

    Args:
        dnase_obj (class): pytest obj
    """

    start = timeit.default_timer()

    dnase_obj.dnase_finder("1:10180", "1:800000")

    end = timeit.default_timer()

    assert 1 > (end - start)


@pytest.mark.skip(reason="최초 데이터 디렉터리 구성할 때 한번만 실행하면 됨")
def test_parsing_wig(dnase_obj):
    """wig file을 분산해서 저장하는 걸 테스트

    Args:
        dnase_obj (class): pytest obj
    """
    dnase_obj.parsing_wig(dnase_obj.dnase_path)

    files = glob.glob("/data/gflas-knockout-efficiency/data/WIG/*.pk")

    assert 118054 == len(files)  ## wig lines start with '#'


def test_dnasetool(dnase_obj):
    """zero-based calculation
    ex> [0.24601, 0.24601, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    Args:
        dnase_obj (class): pytest obj
    """
    input = ("1:10180", "1:10190")
    exp = ([0.24601] * 2) + ([0] * 8)

    pred = dnase_obj.dnase_finder(input[0], input[1])
    print(pred)
    assert exp == pred
