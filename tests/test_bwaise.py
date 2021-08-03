#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

__all__ = ('test_bwaise',)

def test_bwaise():
    from numpy.testing import assert_allclose
    from exposan import bwaise as bw

    assert_allclose(bw.teaA.NPV, -22732728.213841617, rtol=1e-2)
    assert_allclose(bw.teaB.NPV, -2231960.180530764, rtol=1e-2)
    assert_allclose(bw.teaC.NPV, -94340457.26937833, rtol=1e-2)

    assert_allclose(bw.lcaA.total_impacts['GlobalWarming'], 146386354.78674603, rtol=1e-2)
    assert_allclose(bw.lcaB.total_impacts['GlobalWarming'], 11987512.603749081, rtol=1e-2)
    assert_allclose(bw.lcaC.total_impacts['GlobalWarming'], 56832681.43125323, rtol=1e-2)


if __name__ == '__main__':
    test_bwaise()