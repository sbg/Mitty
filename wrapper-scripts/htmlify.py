import os
import base64

html = "<html><head></head><body>"

with open('filelist.txt') as fp:
  for fn in (ffn.strip() for ffn in fp):
    with open(fn, "rb") as image_file:
      im_type = os.path.splitext(os.path.basename(fn))[-1][1:]
      html += '<p align="center"><img src="data:image/{};base64,{}" /></p>'.format(
        im_type, base64.b64encode(image_file.read()).decode('ascii'))

html += "</body>"

with open('report.html', 'w') as fp:
  fp.write(html)