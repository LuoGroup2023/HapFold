#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include <bits/stdc++.h>

KSEQ_INIT(gzFile, gzread)

typedef struct Porec
{
  char *seq_name;
  char *seq_s;
  char *qual;
} Porec;

typedef struct Bam
{
  int start;
  int end;
} Bam;

std::map<std::string, Porec> porec;
std::map<std::string, std::vector<Bam>> bam_all;

int main(int argc, char *argv[])
{
  gzFile fp;
  kseq_t *seq;
  fp = gzopen(argv[1], "r");
  FILE *file_ptr = fopen(argv[2], "r");
  int l;
  if (argc != 3)
  {
    fprintf(stderr, "Usage: %s <in.fasta> <readsinfo.txt>\n", argv[0]);
    return 1;
  }
  long int num = 0;
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0)
  {
    Porec porec_entry;
    porec_entry.seq_name = strdup(seq->name.s);
    porec_entry.seq_s = strdup(seq->seq.s);
    porec_entry.qual = strdup(seq->qual.s);
    std::string seq_id(porec_entry.seq_name);
    size_t dotPosition = seq_id.find(".");

    porec[porec_entry.seq_name] = porec_entry;

    num++;
    if (num % 10000 == 0)
      std::cout << "Loaded " << num << " sequences" << std::endl;
    // porec[]
  }

  char buffer[100];
  if (file_ptr == NULL)
  {
    perror("Error opening file");
    return 1;
  }
  num = 0;
  while (fgets(buffer, sizeof(buffer), file_ptr) != NULL)
  {
    char prefix[30];
    int num1, num2, num3;

    if (sscanf(buffer, "%[^:]:%d:%d", prefix, &num2, &num3) == 3) // Modify the format to adapt to the new input
    {
      printf("Prefix: %s, Start: %d, End: %d\n", prefix, num2, num3);

      Bam bam_p;
      bam_p.start = num2;   // Start position
      bam_p.end = num3 - 1; // End position (assuming 'end' is inclusive, hence subtract 1)

      // Check if the sequence ID exists in the porec map
      if (porec.find(prefix) == porec.end())
      {
        // std::cout << "None in now Porec" << std::endl;
        continue; // Skip if the sequence ID does not exist
      }

      bam_all[prefix].push_back(bam_p); // Append bam_p to bam_all

      num++;
      if (num % 10000 == 0)
        std::cout << "Loaded " << num << " alignment records" << std::endl; // Print progress
    }

    else
    {
      printf("Failed to parse line: %s\n", buffer);
    }
  }

  std::ofstream hic1File("/home/work/yichen/work/2025-pstools/chicken/work/porec_data/porec2hic1.fastq");
  std::ofstream hic2File("/home/work/yichen/work/2025-pstools/chicken/work/porec_data/porec2hic2.fastq");
  num = 0;
  long int size = 0;
  for (const auto &pair : bam_all)
  {
    std::vector<std::string> seq_1;
    std::vector<std::string> qual_1;
    std::string index = pair.first;
    // std::cout << index << std::endl;
    if (num % 10000 == 0)
    {
      std::cout << "Processed " << num << " sequence records" << std::endl;
    }
    // if (index >= 0 && index < static_cast<int>(porec.size())) // Check validity
    //{
    const std::vector<Bam> &bamList = pair.second;
    // std::cout << bamList.size() << std::endl;
    for (const Bam &bamItem : bamList)
    {
      int start = bamItem.start;
      int end = bamItem.end;
      int length = end - start;
      // std::cout << "  Start: " << start << ", End: "<<end <<" index "<<index<< std::endl;
      std::string seq_same(porec[index].seq_s + start, length);
      // std::cout << seq_same << std::endl;
      seq_1.push_back(seq_same);
      std::string seq_qual(porec[index].qual + start, length);
      // std::cout << seq_qual << std::endl;
      qual_1.push_back(seq_qual);
    }

    if (!hic1File.is_open() || !hic2File.is_open())
    {
      std::cerr << "Failed to open one or both files." << std::endl;
      return 1;
    }
    // std::cout << seq_1.size() << std::endl;
    if (seq_1.size() < 2)
      continue;
    int numb = 1;
    size += seq_1.size();
    for (size_t i = 0; i < seq_1.size(); ++i)
    {
      for (size_t j = i + 1; j < seq_1.size(); ++j)
      {
      
        if (i != j)
        {

          hic1File << "@" << porec[index].seq_name << "_" << numb << " " << 1 << std::endl;
          hic1File << seq_1[i] << std::endl;
          hic1File << "+" << std::endl;
          hic1File << qual_1[i] << std::endl;
          // Write to hic2 file
          hic2File << "@" << porec[index].seq_name << "_" << numb << " " << 2 << std::endl;
          hic2File << seq_1[j] << std::endl;
          hic2File << "+" << std::endl;
          hic2File << qual_1[j] << std::endl;
          numb++;
        }
      }
    }
    if (num % 1000 == 0)
    {
      std::cout << "Total alignments generated: " << size << std::endl;
      std::cout << "........................." << std::endl;
    }
    num++;
  }

  hic1File.close();
  hic2File.close();

  fclose(file_ptr);
  kseq_destroy(seq);
  gzclose(fp);
  return 0;
}