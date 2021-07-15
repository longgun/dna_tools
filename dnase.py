import glob
import os
import pickle
from collections import defaultdict


class DnaseTools:
    def __init__(self, dnase_path=""):
        """최소한의 경로만 가진다.
            wig 파일의 경로만 입력으로 받아 처리함

        Args:
            dnase_path (str, optional): wig file path. Defaults to "".
        """
        if not dnase_path:
            dnase_path = "/data/gflas-knockout-efficiency/data/ENCFF185KBK.wig"

        self.dnase_path = dnase_path

        self.file_path = (
            "/data/gflas-knockout-efficiency/data/WIG/dnase_hek293t_"
        )

    def parsing_wig(self, path: str):
        """wig 파일 파싱-저장과 불러오기를 구분하기 위해서 만든 편의상의 메서드

        Args:
            path (str): wig file path
        """
        self.get_data(path)

    def get_file_name(self, line: str) -> str:
        """wig 파일의 #으로 시작하는 라인을 파싱해서 파일이름에 저장되는 구간을 포함시키는 메서드

        Args:
            line (str): #으로 시작하는 라인
            ex> #bedGraph section chr1:0-723230

        Returns:
            str: 쪼개서 저장할 파일들의 이름
            ex> dnase_hek293t_X_98433075_98586930.pk
        """
        res = ""

        section = line.strip().split(" ")[-1]
        chrom, pos = section.split(":")
        chrom = chrom.replace("chr", "")

        start_pos, end_pos = pos.split("-")
        start_pos = int(start_pos)
        end_pos = int(end_pos)

        res = f"{self.file_path}{chrom}_{start_pos}_{end_pos}.pk"
        return res

    def save_res(self, res: dict, file_name: str):
        """dictionary를 pickle로 저장하는 메서드

        Args:
            res (dict): 구간으로 나뉜 wig lines
            file_name (str): 구간정보를 포함한 저장할 파일 이름
        """
        with open(file_name, "wb") as file_handle:
            pickle.dump(res, file_handle)

    def get_data(self, dnase_path: str):
        """실제로 파싱을 수행함 parsing main
            파싱 중간중간마다 res에 저장된 딕셔너리 객체를 저장하고 비움

        Args:
            dnase_path (str): dnase wig file path
        """
        res = defaultdict(dict)
        file_name = ""
        with open(dnase_path) as dnase_fh:
            col = ["chrom", "start", "end", "fold"]
            for line in dnase_fh:
                if line[0] == "#":
                    if res:
                        self.save_res(res, file_name)
                        res = defaultdict(dict)
                    file_name = self.get_file_name(line)
                    continue

                row = dict(zip(col, line.strip().split("\t")))

                pos = f"{row['start']}-{row['end']}"
                res[row["chrom"]][pos] = row["fold"]

        self.save_res(res, file_name)  # 마지막 저장

    def load_data(self, chrom: str, start: int, end: int) -> dict:
        """인풋으로 받은 구간의 데이터를 로딩받아 하나의 딕셔너리로 반환함
        인풋 구간의 데이터가 여러 피클 객체로 나눠져 있을 수 있으므로 딕셔너리를 합치는 부분이 관건

        Args:
            chrom (str): chrN의 형태
            start (int): 0-based start position
            end (int): 0-based end position

        Returns:
            dict: start - end 구간을 모두 포함하는 딕셔너리 객체
        """
        res = dict()

        files = glob.glob(f"{self.file_path}*")

        filename_col = ["dnase", "cellline", "chr", "start", "end"]
        for filename in files:
            splitted = filename.split(".")[0].split("_")
            if len(splitted) != len(filename_col):
                continue
            file_data = dict(zip(filename_col, splitted))
            if file_data["chr"] in chrom:
                if int(file_data["start"]) <= start <= int(file_data["end"]):
                    with open(filename, "rb") as filehandle:
                        res.update(pickle.load(filehandle))
                if int(file_data["start"]) <= end <= int(file_data["end"]):
                    with open(filename, "rb") as filehandle:
                        res.update(pickle.load(filehandle))

        return res

    def dnase_finder(self, pos_from: str, pos_to: str) -> list:
        """문자열 형태의 쿼리를 두개 받아서 dnase value를 list형태로 반환

        Args:
            pos_from (str): chr:pos (0-based)
            pos_to (str): chr:pos (0-based)

        Returns:
            [list]: 주어진 구간의 dnase value pos-by-pos
        """
        res = list()

        start_chrom, start = pos_from.split(":")
        end_chrom, end = pos_to.split(":")

        assert start_chrom == end_chrom, "only work in one chromosome!"

        chrom = f"chr{start_chrom}"
        start = int(start)
        end = int(end)

        assert start < end, "start must be small than end"

        data = self.load_data(chrom, start, end)
        data = data[chrom]

        for key, value in data.items():
            key_start, key_end = key.split("-")

            key_start = int(key_start)
            key_end = int(key_end)

            value = float(value)

            if key_start <= start < key_end:
                for i in range(key_end - start):
                    if i < end - start:
                        res.append(value)

            if key_start <= end <= key_end:
                for i in range(end - key_start):
                    res.append(value)

        return res
