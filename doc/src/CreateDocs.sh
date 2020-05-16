#!/bin/bash

for ((index=0; index<2; index++))
do
    pdflatex manual.tex
    pdflatex IntegralStucture.tex
    pdflatex CCDIS.tex
    pdflatex Transversity.tex
    pdflatex TMDEvolution.tex
    pdflatex DrellYanTMD.tex
    pdflatex SIDISTMD.tex
    pdflatex CovarianceMatrix.tex
    pdflatex DIS.tex
    pdflatex GPDs.tex
    pdflatex MatchingConditions.tex
    pdflatex MSbarMass.tex
done

mv *.pdf ../pdf/.
rm *.aux *.log *.out *.toc *.idx
