from pathlib import Path
from datetime import datetime, timezone
import hashlib, zipfile, xml.etree.ElementTree as ET
root=Path(__file__).resolve().parent.parent
release_date=datetime.now(timezone.utc).strftime('%Y%m%d')
name=f'Lerager-Scripts-{release_date}-v287.zip'
package=root/'updates'/name
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED) as z:
    for file in sorted([*(root/'src').rglob('*.js'), *(root/'doc').rglob('*.html')]):
        z.write(file,file.relative_to(root).as_posix())
sha1=hashlib.sha1(package.read_bytes()).hexdigest()
xri=ET.Element('xri',{'version':'1.0'})
desc=ET.SubElement(xri,'description'); ET.SubElement(desc,'p').text='PixInsight scripts by Morten Lerager: AstroClean and AlignRGB Preview.'
platform=ET.SubElement(xri,'platform',{'os':'all','arch':'noarch','version':'1.9.6:2.99.99'})
pkg=ET.SubElement(platform,'package',{'fileName':name,'sha1':sha1,'type':'script','releaseDate':release_date})
ET.SubElement(pkg,'title').text='Lerager Scripts - AstroClean 2.8.7 / AlignRGB Preview 0.3'
pdesc=ET.SubElement(pkg,'description')
ET.SubElement(pdesc,'p').text='Installs AstroClean and AlignRGB Preview under Script > Lerager.'
ET.indent(xri,space='   ')
(root/"updates"/"updates.xri").write_bytes(b'<?xml version="1.0" encoding="UTF-8"?>\n'+ET.tostring(xri,encoding="utf-8"))
print('NOTE: Repository is unsigned. Signing requires a registered PixInsight signing identity.')
