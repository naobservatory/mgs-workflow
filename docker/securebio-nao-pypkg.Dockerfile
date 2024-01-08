FROM pandas/pandas:pip-all

RUN pip install pandas
RUN pip install biopython
RUN pip install ncbi-genome-download
RUN python -m pip install ete3 six ncbi-genome-download
