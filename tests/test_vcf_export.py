import tempfile
import unittest
from pathlib import Path

import pysam
from seq import io


class TestVcfExport(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.fai = self.root / 'reference.fai'
        self.fai.write_text('chr1\t10000\t0\t50\t51\n')

    def test_header_can_be_parsed_before_adding_sample(self):
        path = self.root / 'template.vcf'
        template = io.get_vcf_format_variant_file(str(path), [io.Chr('chr1', 10000)])
        self.addCleanup(template.close)
        self.assertEqual(list(template.header.samples), ['SAMPLE'])
        self.assertIn('GT', template.header.formats)

    def test_breakpoint_roundtrip(self):
        for svtype in ('DEL', 'DUP', 'INV'):
            for start, end in ((0, 100), (100, 200), (9800, 9999)):
                with self.subTest(svtype=svtype, start=start, end=end):
                    bed = self.root / 'calls.bedpe'
                    vcf = self.root / 'calls.vcf'
                    bed.write_text('chr1\t{}\t{}\tchr1\t{}\t{}\t{}\t0.95\t0/1\n'.format(
                        start, start + 1, end, end + 1, svtype))
                    io.bed2vcf(str(bed), str(vcf), str(self.fai))
                    with pysam.VariantFile(str(vcf)) as output:
                        record = next(output)
                        self.assertEqual(record.pos, start + 1)
                        self.assertEqual(record.stop, end + 1)
                        self.assertEqual(record.info['SVLEN'], end - start)
                        self.assertEqual(record.samples['SAMPLE']['GT'], (0, 1))
                    roundtrip = list(io.vcf_iter(str(vcf)))
                    self.assertEqual(len(roundtrip), 1)
                    self.assertEqual(roundtrip[0].intervalA.start, start)
                    self.assertEqual(roundtrip[0].intervalB.start, end)
                    self.assertFalse(Path(str(vcf) + '.format').exists())

    def test_empty_calls_still_write_valid_vcf(self):
        bed = self.root / 'empty.bedpe'
        vcf = self.root / 'empty.vcf'
        bed.write_text(io.BedRecord.get_bedpe_aux_header() + '\n')
        io.bed2vcf(str(bed), str(vcf), str(self.fai))
        with pysam.VariantFile(str(vcf)) as output:
            self.assertEqual(list(output), [])
            self.assertEqual(list(output.header.samples), ['SAMPLE'])


if __name__ == '__main__':
    unittest.main()
