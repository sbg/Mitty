import time
import pysam
import logging
import sys
from multiprocessing import Process, Queue
from mitty.simulation.sequencing.writefastq import writer, load_qname_sidecar, parse_qname
from mitty.simulation.readgenerate import DNA_complement
logger = logging.getLogger(__name__)
__process_stop_code__ = 'SETECASTRONOMY'
in_queue, out_queue = Queue(10000), Queue(10000)
processed_tuples={}


def read2vlist(read):
  vlist=""
  for cigar in read.cigartuples:
    if cigar[0]==1:
     vlist=vlist+","+str(cigar[1])
    elif cigar[0]==2:
      vlist=vlist+","+str(-cigar[1])
    elif cigar[0]==8:
      vlist=vlist+",0"
  if len(vlist)>1:
    vlist="("+vlist[1:]+")"
  return vlist

def bam_to_truth(bam_fname_in,mq_threshold, sample_name,output_prefix):
  logger.debug('Starting filtering ...')
  t0 = time.time()
  wr = Process(target=writer, args=(output_prefix+"-1.fq", output_prefix+"-lq.txt", output_prefix+"-2.fq", out_queue))
  wr.start()
  i=0
  for read in pysam.AlignmentFile(bam_fname_in, "rb"):  
    i=i+1
    if i % 32768==0:
      logger.debug('Processed {} reads'.format(i))
    if read.is_paired:
      if read.qname in processed_tuples:
        #mate 2
        if (read.mapping_quality >= mq_threshold) and (read.next_reference_id==read.reference_id) and (read.reference_name==processed_tuples[read.qname][2]): 
          if read.is_reverse:
            out_queue.put((None,sample_name,read.reference_name,0,((0 ,read.pos+1,read.cigarstring,read2vlist(read),"",read.seq.translate(DNA_complement)[::-1],read.qual[::-1]),(0 ,processed_tuples[read.qname][4]+1,processed_tuples[read.qname][5],processed_tuples[read.qname][6],"",processed_tuples[read.qname][1],processed_tuples[read.qname][0]))))
          else:
            out_queue.put((None,sample_name,read.reference_name,0,((0 ,read.pos+1,read.cigarstring,read2vlist(read),"",read.seq,read.qual),(0 ,processed_tuples[read.qname][4]+1,processed_tuples[read.qname][5],processed_tuples[read.qname][6],"",processed_tuples[read.qname][1],processed_tuples[read.qname][0]))))
        del processed_tuples[read.qname]

      else:
        #mate 1
        if (not (read.is_unmapped or read.mate_is_unmapped) and (read.mapping_quality >= mq_threshold)) and (read.next_reference_id==read.reference_id): 
          if read.is_reverse:
            processed_tuples[read.qname]=(read.qual[::-1],read.seq.translate(DNA_complement)[::-1],read.reference_name,read.is_reverse,read.pos,read.cigarstring,read2vlist(read))
          else:
            processed_tuples[read.qname]=(read.qual,read.seq,read.reference_name,read.is_reverse,read.pos,read.cigarstring,read2vlist(read))
          
    else:
      #not pair-end
      if not read.is_unmapped and (read.mapping_quality >= mq_threshold) and (read.next_reference_id==read.reference_id):
        if read.is_reverse:
          out_queue.put((None,sample_name,read.reference_name,0,(( 0 ,read.pos+1,read.cigarstring,read2vlist(read),"",read.seq.translate(DNA_complement)[::-1],read.qual[::-1]),)))
        else:
          out_queue.put((None,sample_name,read.reference_name,0,(( 0 ,read.pos+1,read.cigarstring,read2vlist(read),"",read.seq,read.qual),)))
          
  out_queue.put(__process_stop_code__)
  wr.join()
  t1 = time.time()
  logger.debug('{} reads were processed'.format(i))
  logger.debug('Took {} s'.format(t1 - t0))
