#!/bin/bash

CPPFLAGS="-I${PREFIX}/include"
LDFLAGS="-L${PREFIX}/lib"

echo "### Create dekpul conda folder ###"
mkdir -p "$PREFIX"/share/dekupl/bin
mkdir -p "$PREFIX"/share/dekupl/share

echo "### Install computeNF ###"
rm -f share/computeNF/computeNF && make -C share/computeNF LIBS="-lz -lm $CPPFLAGS $LDFLAGS"
cp share/computeNF/computeNF "$PREFIX"/share/dekupl/bin

echo "### Install joinCounts ###"
rm -f share/joinCounts/joinCounts && make -C share/joinCounts LIBS="-lz -lm $CPPFLAGS $LDFLAGS"
cp share/joinCounts/joinCounts "$PREFIX"/share/dekupl/bin

echo "### Install mergeTags ###"
rm -f share/mergeTags/mergeTags && make -C share/mergeTags LIBS="-lz -lm $CPPFLAGS $LDFLAGS"
cp share/mergeTags/mergeTags "$PREFIX"/share/dekupl/bin

echo "### Install TtestFilter ###"
rm -f share/TtestFilter/TtestFilter && make -C share/TtestFilter LIBS="-lz -lm $CPPFLAGS $LDFLAGS"
cp share/TtestFilter/TtestFilter "$PREFIX"/share/dekupl/bin

echo "### Install Kalisto ###"
wget https://github.com/pachterlab/kallisto/releases/download/v0.43.0/kallisto_linux-v0.43.0.tar.gz -O kallisto.tar.gz
tar -xzf kallisto.tar.gz
cp kallisto_linux-v0.43.0/kallisto "$PREFIX"/share/dekupl/bin

echo "### Install dekupl ###"
cp bin/DESeq2_diff_method.R "$PREFIX"/share/dekupl/bin
cp bin/DESeq2_ref_transcripts.R "$PREFIX"/share/dekupl/bin
cp bin/diffFilter.pl "$PREFIX"/share/dekupl/bin
cp bin/mergeCounts.pl "$PREFIX"/share/dekupl/bin
cp bin/revCompFastq.pl "$PREFIX"/share/dekupl/bin
cp bin/Ttest_diff_method.R "$PREFIX"/share/dekupl/bin

cp Snakefile "$PREFIX"/share/dekupl
cp config.json "$PREFIX"/share/dekupl

cp dekupl-run.sh "$PREFIX"/share/dekupl

echo "### Install bin ###"
ln -svf "$PREFIX"/share/dekupl/dekupl-run.sh "$PREFIX"/bin/dekupl-run

