from bs4 import BeautifulSoup as BS
import urllib3
import time

url = "http://origin.tubic.org/deg/public/index.php/organism/eukaryotes/DEG2011.html?lineage=eukaryotes&id=DEG2011&page="

with open("essGene.csv", "a+") as f:
    f.write("geneName,DEG_id,uniprot_id\n")

for i in range(1, 55):
    http_pool = urllib3.connection_from_url(url+str(i))
    r = http_pool.urlopen('GET',url+str(i))
    html = r.data.decode('utf-8')
    soup = BS(html, features="lxml")
    table = soup.find(id='browse').select(".table")[0]
    # loop each protein row
    for item in table.find_all('tr'):
        if str(item.parent.name) == "thead": continue

        # get gene name
        cols = item.find_all("td")
        geneName = cols[2].text
        uniprot_id = None

        # find uniprot id
        DEG_id = cols[1].find_all('a')[0].text
        id_url = "http://origin.tubic.org/deg/public/index.php/information/eukaryotes/{}.html".format(DEG_id)
        id_r = urllib3.connection_from_url(id_url).urlopen('GET', id_url)
        id_html = id_r.data.decode('utf-8')
        id_soup = BS(id_html, features="lxml")

        # get uniprot id
        idTable = id_soup.select(".table")[0]
        for idRow in idTable.find_all('tr'):
            if "Uniprot" not in idRow.text: continue
            uniprot_id = idRow.find_all('td')[1].text

        print(geneName, DEG_id, uniprot_id)
        with open("essGene.csv", "a+") as f:
            f.write("{},{},{}\n".format(geneName, DEG_id, uniprot_id))
        time.sleep(1)

    time.sleep(1)