#! /bin/sh
echo "###################################################"
echo "Checking integrate --interactive-integrands..."
echo "###################################################"
./test-interactive-integrands.pl "--cone-decompose" || exit 1
./test-interactive-integrands.pl "--triangulate" || exit 1
