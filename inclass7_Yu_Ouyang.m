%Inclass assignment 7. 
% 1. The gene Cdx2 is a crucial transcription factor involved in number of
% developmental stages. Use the UCSC genome browser to answer the following questions
% about it:

% A. What human chromosome is it located on?
% Yu Ouyang's answer: It is on chromosome 13

% B. How many exons does it have?
% Yu Ouyang: 3 exons

% C. What is the precise position of its stop codon in the genome?
% Yu Ouyang: 27,969,006

% D. Identify at least one difference in sequence between human and mouse
% CDX2.
% Yu Ouyang: Size difference: Human gene Cdx2 is 5892bp for coding region
% and 9003bp for transcripts with UTRs, while mouse gene Cdx2 is 5141bp for
% coding region and 6350bp for transcripts with UTRs. 

% E. In which human tissues is it expressed most abundantly?
% Yu Ouyang: Colon-Transverse

%2. A. Use the unigene database to find the accession number for a genbank
% entry containing the complete coding sequence of Cdx2. 
% Yu Ouyang: NM_001265

% B. Use MATLAB to read the genbank information corresponding to that
% accession number. 
fasta_dat = fastaread('sequence.fasta.txt');
% or
genbank_dat = genbankread('sequence.gb.txt');
accession = genbank_dat.Accession;
genbank_dat2 = getgenbank(accession);


% C. Use the information read in to find the position of the start and stop
% codon within the sequence. What are the parts of the sequence before the start codon 
% and after the stop codon?
start_pos = strfind(fasta_dat.Sequence, 'GCTCCCGGACCCTCGCCACCATG')+20;
stop_codon = [strfind(fasta_dat.Sequence, 'TAA') strfind(fasta_dat.Sequence, 'TGA') strfind(fasta_dat.Sequence, 'TAG')];
Olengths = stop_codon - start_pos+3;
good_length = 1e8;
good_ind = 0;
for jj = 1:length(Olengths)
    if Olengths(jj) > 0 && ...
            mod(Olengths(jj),3) == 0 && ...
            Olengths(jj) < good_length
        good_length = Olengths(jj);
        good_ind = jj;
    end
end
stop_pos = stop_codon(good_ind);
strfind(fasta_dat.Sequence, 'CAAAAAAA');
% Yu Ouyang: The start codon is at 363, and the stop codon is at 1302. The
% 'CAAAAAA' poly-A tail is at 2351. In the structure CDS, it says the start and
% stop codon positions are 363 and 1304. 

% D. Use the protein_id to read the information on the protein. Use the
% information you read to determine where the homeobox domain of the protein is.
% Hint: see the field "Features". 

protein_dat = getgenpept(genbank_dat2.CDS.protein_id);
info = protein_dat.Features;
% Yu Ouyang: The "Homeobox" is at 190-242. In "Region", 190...242, it says
% '/region_name = "Homeobox"'.