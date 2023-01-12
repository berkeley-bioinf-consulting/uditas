library(ggplot2)
library(readxl)
library(officer)

add_to_ppt = function(ppt, plt, title) {
	    ppt = add_slide(ppt)
    ppt = ph_with(x=ppt, value=title, location = ph_location_type(type = "title"))
        ppt = ph_with(x=ppt, value=plt, location = ph_location_type(type = "body"), bg = "transparent")
        ppt
}

# open big results file
df = readxl::read_excel('all_results/data_big_results.xlsx')
df[df$Sample == 'S8_maybe', 'Sample'] = 'S8'
# plot read counts, collapsed read counts, fraction on-target
ppt = read_pptx()

gpt = ggplot(df, aes(Sample, read_count)) + geom_col() + ylab('Raw Read Count')
ppt = add_to_ppt(ppt, gpt, 'No obvious library bias in read count')
gpt = ggplot(df, aes(Sample, percent_aligned)) + geom_col() + ylab('% Aligned')
ppt = add_to_ppt(ppt, gpt, 'High variability in mapped read percentage')
gpt = ggplot(df, aes(Sample, total_aligned_collapsed)) + geom_col() + ylab('Unique Aligned Read Count')
ppt = add_to_ppt(ppt, gpt, 'Read counts after UMI deduplication are very consistent – even split sample')
gpt = ggplot(df, aes(Sample, total_aligned_amplicons_collapsed)) + geom_col() + ylab('Unique On-target Read Count')
ppt = add_to_ppt(ppt, gpt, 'Unique on-target read counts after UMI deduplication are very variable')
gpt = ggplot(df, aes(Sample, total_aligned_amplicons_collapsed / total_aligned_collapsed * 100)) + geom_col() + ylab('% Unique On-target')
ppt = add_to_ppt(ppt, gpt, 'On-target percentage of unique reads')
gpt = ggplot(df, aes(Sample, total_aligned_amplicons_collapsed, fill=primer)) + geom_col() + ylab('Unique On-target Read Count')
ppt = add_to_ppt(ppt, gpt, 'On-target read counts don’t seem to by driven by which primer was used')
gpt = ggplot(df, aes(Sample, percent_aligned_all_amplicons)) + geom_col() + ylab('% Raw On-target')
ppt = add_to_ppt(ppt, gpt, 'Raw on-target percentage (as reported by Uditas) is probably not meaningful)')
print(ppt, target='uditas_QC.pptx')

