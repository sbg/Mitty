import os
import time
import logging


logger = logging.getLogger(__name__)


def read_fastq(fp, ipe, f_size=None, max_templates=None):
  """

  :param fp: pointer to file/stream
  :param ipe: True if this is an interleaved paired end file
  :param f_size:  os.stat(fastq).st_size. Leave as None if unknowable, like a stream
  :param max_templates: quit after these many templates read
  :return:
  """
  template_notification_interval = 100000
  template_counter = 0
  t0 = time.time()

  template = []
  read = []
  reads_in_template_count = 2 if ipe else 1
  line_cnt = 4
  for ln in fp:
    read.append(ln[:-1])
    line_cnt -= 1
    if line_cnt == 0:
      template.append(read)
      read = []
      line_cnt = 4
      reads_in_template_count -= 1
      if reads_in_template_count == 0:
        # in_queue.put(template)
        yield template
        template_counter += 1
        if template_counter % template_notification_interval == 0:
          progress_message(template_counter, t0, fp, f_size)
        if max_templates is not None and template_counter >= max_templates:
          logger.debug('Stopping at {} templates, as asked for'.format(max_templates))
          break
        template = []
        reads_in_template_count = 2 if ipe else 1

  t1 = time.time()
  logger.debug('Took {:0.5}s to process {} templates ({:0.7} tmpl/s)'.format(t1 - t0, template_counter, template_counter/(t1 - t0)))


def progress_message(template_counter, t0, fp, f_size=None):
  progr = ', {:0.5}/{:0.5} MB'.format(fp.tell() / 1e6, f_size /1e6) if f_size is not None else ''
  logger.debug(
    'Processed {} templates '
    '({:0.7} templates/s{})'.format(template_counter, template_counter / (time.time() - t0), progr))
