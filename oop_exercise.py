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

    print(a.describe())  
    print(a.length())  
    print(a.overlaps(b))  
    print(a.overlaps(c))  

#TASK 2
class Exon(GenomicFeature):
    def __init__(self, chromosome, start, end, strand, exon_number):
        super().__init__(chromosome, start, end, strand)
        if not isinstance(exon_number, int):
            raise ValueError("Exon number is incorrect")

        self.exon_number = exon_number

    def describe(self) -> str:
        return f"{type(self).__name__} {self.chromosome}:{self.start}-{self.end}({self.strand}) exon #{self.exon_number}"

if __name__ == "__main__":
    features = [
        GenomicFeature("chr1", 1000, 5000, "+"),
        Exon("chr1", 1000, 1200, "+", 1),
        Exon("chr1", 3000, 3300, "+", 2),
    ]
    for feature in features:
        print(feature.describe())

#TASK 3
class Exon(GenomicFeature):
    def __init__(self, chromosome, start, end, strand, exon_number):
        super().__init__(chromosome, start, end, strand)
        if not isinstance(exon_number, int):
            raise ValueError("Exon number is incorrect")

        self.exon_number = exon_number

    def describe(self) -> str:
        return f"{type(self).__name__} {self.chromosome}:{self.start}-{self.end}({self.strand}) exon #{self.exon_number}"

class Gene(GenomicFeature):
    def __init__(self, chromosome, start, end, strand, name):
        super().__init__(chromosome, start, end, strand)
        if not isinstance(name, str):
            raise ValueError("Name is incorrect")
        self.name = name
        self.exons = []

    def add_exon(self, exon):
        if not isinstance(exon, Exon):
            raise ValueError ("Exon is incorrect")
        self.exons.append(exon)
    def total_exon_length(self) -> int:
        return sum(exon.length() for exon in self.exons)
    def describe(self) -> str:
        return f"{type(self).__name__} {self.chromosome}:{self.start}-{self.end}({self.strand}) {self.name}, {len(self.exons)} exon(s)"

class Variant(GenomicFeature):
    def __init__(self, chromosome, start, end, strand, ref_allele, alt_allele):
        super().__init__(chromosome, start, end, strand)
        if not isinstance(ref_allele, str):
            raise ValueError("Ref allele is incorrect")
        if not isinstance(alt_allele, str):
            raise ValueError("Alt allele is incorrect")

        self.ref_allele = ref_allele
        self.alt_allele = alt_allele

    def variant_type(self) -> str:
        if len(self.ref_allele) ==1 and len(self.alt_allele) == 1:
            return "SNP"
        if len(self.ref_allele) < len(self.alt_allele):
            return "insertion"
        if len(self.ref_allele) > len(self.alt_allele):
            return "deletion"
        else:
            return "MNV"

    def describe(self) -> str:
        return f"{type(self).__name__} {self.chromosome}:{self.start}-{self.end}({self.strand}) {self.ref_allele}>{self.alt_allele} ({self.variant_type()})"
if __name__ == "__main__":
    genes: dict[str, Gene] = {}
    variants: list[Variant] = []

    with open("oop_data.tsv") as file:
        for line in file:
            line = line.strip()

            if not line or line.startswith("type") or line.startswith("#"):
                continue

            parts = line.split("\t")
            feature_type = parts[0]
            chromosome = parts[1]
            start = int(parts[2])
            end = int(parts[3])
            strand = parts[4]
            field_a = parts[5]
            if len(parts) > 6:
                field_b = parts[6]
            else:
                field_b = ""

            if feature_type == "gene":
                genes[field_a] = Gene(chromosome, start, end, strand, field_a)
            elif feature_type == "exon":
                exon = Exon(chromosome, start, end, strand, int(field_b))
                if field_a in genes:
                    genes[field_a].add_exon(exon)
            elif feature_type == "variant":
                variants.append(
                    Variant(chromosome, start, end, strand, field_a, field_b)
                )

    all_features = list(genes.values()) + variants
    for feature in all_features:
        if isinstance(feature, Gene):
            print(
                f"{feature.describe()} | Total exon length: {feature.total_exon_length()}"
            )
        else:
            print(feature.describe())

    for var in variants:
        found_in_gene = False
        for gene in genes.values():
            if var.overlaps(gene):
                found_in_gene = True
                in_exon = any(var.overlaps(exon) for exon in gene.exons)
                if in_exon:
                    print(f"{var.describe()} -> inside {gene.name} (in exon)")
                else:
                    print(
                        f"{var.describe()} -> inside {gene.name} (not in exon)"
                    )
        if not found_in_gene:
            print(f"{var.describe()} -> intergenic")