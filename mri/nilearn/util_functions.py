def createSubjectlist(subjects):

	subjList = []
	
	for subj in subjects:
		subjList.append('sub-' + str(subj).zfill(2))
	
	return subjList


def create_html_report(image_list, output_filename, title=''):
	
	import base64
	"""
	Currently does not work in python3
	"""
	
	f = open(output_filename, 'w')
	f.write('<html>\n')
	f.write('<h1 style="font-family:helvetica;">')
	if title != '':
		f.write('<head><title>' + title + '</title></head>\n')
	f.write('<body><p><font size="14">' + title + '</font></p></body>\n')
	if title != '':
		f.write('<hr>')
	
	for img in image_list:
		image = base64.b64encode(open(img, 'rb').read()).decode().replace('\n', '')
		f.write('<img src="data:image/png;base64,{0}" width="1650">'.format(image))
		f.write('<hr>')
		f.write('\n')
		
	f.write('</h1>')
	f.write('</html>\n')
	f.close()
	