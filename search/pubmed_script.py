from Bio import Entrez
import time

Entrez.email = "araz.rawshani@gu.se"  # Replace with your email

def fetch_batch(id_list):
    handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    return records

def fetch_articles(term, total_articles, mindate, maxdate):
    data = []
    batch_size = 100  # You can adjust this as needed
    batches = total_articles // batch_size

    for batch in range(batches):
        handle = Entrez.esearch(db="pubmed",
                                retmax=batch_size,
                                retstart=batch * batch_size,
                                term=term,
                                datetype="pdat",
                                mindate=mindate,
                                maxdate=maxdate)
        record = Entrez.read(handle)
        handle.close()

        id_list = record.get("IdList", [])
        if not id_list:
            print(f"No articles found for term '{term}' between {mindate} and {maxdate}.")
            break  # Exit the loop if no IDs were found
        
        records = fetch_batch(id_list)

        for record in records['PubmedArticle']:
            article = record['MedlineCitation']['Article']
            title = article.get('ArticleTitle', 'No title available')
            pmid = record['MedlineCitation']['PMID']
            year = article['Journal']['JournalIssue']['PubDate'].get('Year', 'No year available')

            # Extracting DOI
            doi = None
            if 'ELocationID' in article:
                for location in article['ELocationID']:
                    if isinstance(location, str) and location.startswith('10.'):
                        doi = location
                        break

            # Extracting all affiliations for all authors
            all_affiliations = []
            if 'AuthorList' in article:
                for author in article['AuthorList']:
                    if 'AffiliationInfo' in author:
                        author_affiliations = [info['Affiliation'] for info in author['AffiliationInfo']]
                        all_affiliations.extend(author_affiliations)

            # Combine all unique affiliations for the current article
            unique_affiliations = list(set(all_affiliations))

            data.append({
                'Title': title,
                'PMID': pmid,
                'DOI': doi,
                'Year': year,
                'Affiliations': '; '.join(unique_affiliations)  # Separated by semicolon
            })

        # Pause for a short time before the next batch
        time.sleep(2)

    return data
