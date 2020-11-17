def create_html_report(image_list, output_filename, title=''):
	
	f = open(output_filename, 'w')
	f.write('<html>\n')
	f.write('<h1 style="font-family:helvetica;">')
	if title != '':
		f.write('<head><title>' + title + '</title></head>\n')
	f.write('<body><p><font size="14">' + title + '</font></p></body>\n')
	if title != '':
		f.write('<hr>')
	
	for img in image_list:
		image = open(img, 'rb').read().encode('base64').replace('\n', '')
		f.write('<img src="data:image/png;base64,{0}" width="1900">'.format(image))
		f.write('<hr>')
		f.write('\n')
		
	f.write('</h1>')
	f.write('</html>\n')
	f.close()
