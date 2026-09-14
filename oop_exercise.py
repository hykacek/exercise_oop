#TASK 1

class GenomicFeature:
    def __init__ (self, chromosome, start, end, strand):
        if not isinstance(chromosome, str):
            raise ValueError("Chromosome is incorrect")
        if not isinstance(start, int) and start>=1 and start<=end:
            raise ValueError("Start is incorrect")
        if not isinstance(end, int) and end>=1 and end>=start:
            raise ValueError("End is incorrect")
        if strand not in ("+","-"):
            raise ValueError("Strand is incorrect")

        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.strand = strand

    def length(self) -> int:
        return self.end - self.start + 1
    def overlaps(self, other) -> bool:
        if not isinstance(other, GenomicFeature):
            return False
        if self.chromosome != other.chromosome:
            return False
        return self.start <= other.end and other.start <= self.end

    def describe(self) -> str:
        return f"{type(self).__name__} {self.chromosome}:{self.start}-{self.end}({self.strand})"

if __name__ == "__main__":
    a = GenomicFeature("chr1", 1000, 5000, "+")
    b = GenomicFeature("chr1", 4800, 6000, "+")
    c = GenomicFeature("chr2", 1000, 5000, "+")

    print(a.describe())  # GenomicFeature chr1:1000-5000(+)
    print(a.length())  # 4001
    print(a.overlaps(b))  # True  (4800-5000 shared)
    print(a.overlaps(c))  # False (different chromosome)
    GenomicFeature("chr1", 5000, 1000, "+")
