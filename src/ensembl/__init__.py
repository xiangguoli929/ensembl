from collections import defaultdict
from enum import Enum
from functools import singledispatchmethod
from urllib.parse import urljoin

import requests
from requests.adapters import HTTPAdapter, Retry


class ContentType(Enum):
    json = "application/json"
    xml = "text/xml"
    nh = "text/x-nh"
    phyloxml = "text/x-phyloxml+xml"
    orthoxml = "text/x-orthoxml+xml"
    gff3 = "text/x-gff3"
    fasta = "text/x-fasta"
    bed = "text/x-bed"
    seqxml = "text/x-seqxml+xml"
    text = "text/plain"
    yaml = "text/x-yaml"
    jsonp = "text/javascript"


class Ensembl:
    headers = defaultdict(str)

    def __init__(self, assembly="GRCh38", scheme="http"):
        self.session = requests.Session()
        adapter = HTTPAdapter(max_retries=Retry(backoff_factor=3600/55000,
                                                respect_retry_after_header=True, status_forcelist=[429], allowed_methods=["GET", "POST"]))
        self.session.mount("http://", adapter)
        self.session.mount("https://", adapter)
        match assembly, scheme:
            case "GRCh38", "http":
                self.server = "http://rest.ensembl.org"
            case "GRCh37", "http":
                self.server = "http://grch37.rest.ensembl.org"
            case "GRCh38", "https":
                self.server = "https://rest.ensembl.org"
            case "GRCh37", "https":
                self.server = "https://grch37.rest.ensembl.org"

    def content_type(self, format):
        self.headers["Content-Type"] = ContentType[format].value
        return self.headers

    def get(self, endpoint, params, format):
        response = self.session.get(urljoin(
            self.server, endpoint), headers=self.content_type(format), params=params)
        if response.ok:
            if self.headers["Content-Type"] == "application/json":
                return response.json()
            else:
                return response.text
        else:
            response.raise_for_status()

    def post(self, endpoint, params, json, format):
        response = self.session.post(urljoin(
            self.server, endpoint), headers=self.content_type(format), params=params, json=json)
        if response.ok:
            if self.headers["Content-Type"] == "application/json":
                return response.json()
            else:
                return response.text
        else:
            response.raise_for_status()

    @singledispatchmethod
    def archive(self, id: str, format="json"):
        """Uses the given identifier to return its latest version"""
        return self.get(endpoint=f"archive/id/{id}", params={}, format=format)

    @archive.register
    def _(self, id: list[str], format="json"):
        """Retrieve the latest version for a set of identifiers"""
        return self.post(endpoint=f"archive/id", json={"id": id}, params={}, format=format)

    def cafe_genetree_id(self, id, format="json", compara=None, nh_format=None):
        """Retrieves a cafe tree of the gene tree using the gene tree stable identifier"""
        return self.get(endpoint=f"cafe/genetree/id/{id}", params=dict(compara=compara, nh_format=nh_format), format=format)

    def cafe_genetree_member_id(self, id, format="json", compara=None, nh_format=None, db_type=None, object_type=None, species=None):
        """Retrieves the cafe tree of the gene tree that contains the gene / transcript / translation stable identifier"""
        return self.get(endpoint=f"cafe/genetree/member/id/{id}", params=dict(compara=compara, nh_format=nh_format, db_type=db_type, object_type=object_type, species=species), format=format)

    def cafe_genetree_member_symbol(self, symbol, species="human", format="json", compara=None, db_type=None, external_db=None, nh_format=None, object_type=None):
        """Retrieves the cafe tree of the gene tree that contains the gene identified by a symbol"""
        return self.get(endpoint=f"cafe/genetree/member/symbol/{species}/{symbol}", params=dict(compara=compara, db_type=db_type, external_db=external_db, nh_format=nh_format, object_type=object_type), format=format)

    def genetree_id(self, id, format="json", aligned=None, cigar_line=None, clusterset_id=None, compara=None, nh_format=None, prune_species=None, prune_taxon=None, sequence=None):
        """Retrieves a gene tree for a gene tree stable identifier"""
        return self.get(endpoint=f"genetree/id/{id}", params=dict(aligned=aligned, cigar_line=cigar_line, clusterset_id=clusterset_id, compara=compara, nh_format=nh_format,
                                                                  prune_species=prune_species, prune_taxon=prune_taxon, sequence=sequence), format=format)

    def genetree_member_id(self, id, format="json", aligned=None, cigar_line=None, clusterset_id=None, compara=None, db_type=None, nh_format=None, object_type=None,
                           prune_species=None, prune_taxon=None, sequence=None, species=None):
        """Retrieves the gene tree that contains the gene / transcript / translation stable identifier"""
        return self.get(endpoint=f"genetree/member/id/{id}", params=dict(aligned=aligned, cigar_line=cigar_line, clusterset_id=clusterset_id, compara=compara, db_type=db_type, nh_format=nh_format,
                                                                         object_type=object_type, prune_species=prune_species, prune_taxon=prune_taxon, sequence=sequence, species=species), format=format)

    def genetree_member_symbol(self, species, symbol, format="json", aligned=None, cigar_line=None, clusterset_id=None, compara=None, db_type=None, external_db=None, nh_format=None, object_type=None,
                               prune_species=None, prune_taxon=None, sequence=None):
        """Retrieves the gene tree that contains the gene identified by a symbol"""
        return self.get(endpoint=f"genetree/member/symbol/{species}/{symbol}", params=dict(aligned=aligned, cigar_line=cigar_line, clusterset_id=clusterset_id, compara=compara, db_type=db_type, nh_format=nh_format,
                                                                                           object_type=object_type, external_db=external_db, prune_species=prune_species, prune_taxon=prune_taxon, sequence=sequence), format=format)

    def alignment_region(self, species, region, format="json", aligned=None, compact=None, compara=None, display_species_set=None, mask=None, method=None, species_set=None, species_set_group=None):
        """Retrieves genomic alignments as separate blocks based on a region and species"""
        return self.get(endpoint=f"alignment/region/{species}/{region}", params=dict(aligned=aligned, compact=compact, compara=compara, display_species_set=display_species_set, mask=mask, method=method, species_set=species_set,
                                                                                     species_set_group=species_set_group), format=format)

    def homology_id(self, id, format="json", aligned=None, cigar_line=None, compara=None, sequence=None, target_species=None, target_taxon=None, type=None):
        """Retrieves homology information (orthologs) by Ensembl gene id"""
        return self.get(endpoint=f"homology/id/{id}", params=dict(aligned=aligned, cigar_line=cigar_line, compara=compara, sequence=sequence, target_species=target_species, target_taxon=target_taxon,
                                                                  type=type), format=format)

    def homology_symbol(self, species, symbol, format="json", aligned=None, cigar_line=None, compara=None, external_db=None, sequence=None,
                        target_species=None, target_taxon=None, type=None):
        """Retrieves homology information (orthologs) by symbol"""
        return self.get(f"homology/symbol/{species}/{symbol}", params=dict(aligned=aligned, cigar_line=cigar_line, compara=compara, external_db=external_db,
                                                                           sequence=sequence, target_species=target_species, target_taxon=target_taxon, type=type), format=format)

    def xrefs_symbol(self, species, symbol, format="json", db_type=None, external_db=None, object_type=None):
        """Looks up an external symbol and returns all Ensembl objects linked to it.
        This can be a display name for a gene/transcript/translation, a synonym or an externally linked reference.
        If a gene's transcript is linked to the supplied symbol the service will return both gene and transcript (it supports transient links)."""
        return self.get(f"xrefs/symbol/{species}/{symbol}", params=dict(db_type=db_type, external_db=external_db, object_type=object_type), format=format)

    def xrefs_id(self, id, format="json", all_levels=None, db_type=None, external_db=None, object_type=None, species=None):
        """Perform lookups of Ensembl Identifiers and retrieve their external references in other databases"""
        return self.get(f"xrefs/id/{id}", params=dict(all_levels=all_levels, db_type=db_type, external_db=external_db, object_type=object_type, species=species), format=format)

    def xrefs_name(self, species, name, format="json", db_type=None, external_db=None):
        """Performs a lookup based upon the primary accession or display label of an external reference and returning the information we hold about the entry"""
        return self.get(f"xrefs/name/{species}/{name}", params=dict(db_type=db_type, external_db=external_db), format=format)

    def info_analysis(self, species, format="json"):
        """List the names of analyses involved in generating Ensembl data."""
        return self.get(f"info/analysis/{species}", params={}, format=format)

    def info_assembly(self, species, format="json", bands=None, synonyms=None):
        """List the currently available assemblies for a species, along with toplevel sequences, chromosomes and cytogenetic bands."""
        return self.get(f"info/assembly/{species}", params=dict(bands=bands, synonyms=synonyms), format=format)

    def info_assembly_region_name(self, species, region_name, format="json", bands=None, synonyms=None):
        """Returns information about the specified toplevel sequence region for the given species."""
        return self.get(f"info/assembly/{species}/{region_name}", params=dict(bands=bands, synonyms=synonyms), format=format)

    def info_biotypes(self, species, format="json"):
        """List the functional classifications of gene models that Ensembl associates with a particular species.
        Useful for restricting the type of genes/transcripts retrieved by other endpoints."""
        return self.get(f"info/biotypes/{species}", params={}, format=format)

    def info_biotypes_groups(self, format="json", group=None, object_type=None):
        """Without argument the list of available biotype groups is returned.
        With :group argument provided, list the properties of biotypes within that group.
        Object type (gene or transcript) can be provided for filtering."""
        return self.get(f"info/biotypes/groups", params=dict(group=group, object_type=object_type), format=format)

    def info_biotypes_name(self, name, format="json", object_type=None):
        """List the properties of biotypes with a given name. Object type (gene or transcript) can be provided for filtering."""
        return self.get(f"info/biotypes/name/{name}", params=dict(object_type=object_type), format=format)

    def info_compara_methods(self, format="json", CLASS=None, compara=None):
        """List all compara analyses available (an analysis defines the type of comparative data)."""
        return self.get(f"info/compara/methods", params=dict(CLASS=CLASS, compara=compara), format=format)

    def info_compara_species_sets(self, method, format="json", compara=None):
        """List all collections of species analysed with the specified compara method."""
        return self.get(f"info/compara/species_sets/{method}", params=dict(compara=compara), format=format)

    def info_comparas(self, format="json"):
        """Lists all available comparative genomics databases and their data release."""
        return self.get(f"info/comparas", params={}, format=format)

    def info_data(self, format="json"):
        """Shows the data releases available on this REST server."""
        return self.get(f"info/data", params={}, format=format)

    def info_eg_version(self, format="json"):
        """Returns the Ensembl Genomes version of the databases backing this service"""
        return self.get(f"info/eg_version", params={}, format=format)

    def info_external_dbs(self, species, format="json", feature=None, filter=None):
        """Lists all available external sources for a species."""
        return self.get(f"info/external_dbs/{species}", params=dict(feature=feature, filter=filter), format=format)

    def info_divisions(self, format="json"):
        """Get list of all Ensembl divisions for which information is available"""
        return self.get(f"info/divisions", params={}, format=format)

    def info_genomes(self, name, format="json", expand=None):
        """Find information about a given genome"""
        return self.get(f"info/genomes/{name}", params=dict(expand=expand), format=format)

    def info_genomes_accession(self, accession, format="json", expand=None):
        """Find information about genomes containing a specified INSDC accession"""
        return self.get(f"info/genomes/accession/{accession}", params=dict(expand=expand), format=format)

    def info_genomes_assembly(self, assembly_id, format="json", expand=None):
        """Find information about a genome with a specified assembly"""
        return self.get(f"info/genomes/assembly/{assembly_id}", params=dict(expand=expand), format=format)

    def info_genomes_division(self, division, format="json", expand=None):
        """Find information about all genomes in a given division. May be large for Ensembl Bacteria."""
        return self.get(f"info/genomes/division/{division}", params=dict(expand=expand), format=format)

    def info_genomes_taxonomy(self, taxon_name, format="json", expand=None):
        """Find information about all genomes beneath a given node of the taxonomy"""
        return self.get(f"info/genomes/taxonomy/{taxon_name}", params=dict(expand=expand), format=format)

    def info_ping(self, format="json"):
        """Checks if the service is alive."""
        return self.get(f"info/ping", params={}, format=format)

    def info_rest(self, format="json"):
        """Shows the current version of the Ensembl REST API."""
        return self.get(f"info/rest", params={}, format=format)

    def info_software(self, format="json"):
        """Shows the current version of the Ensembl API used by the REST server."""
        return self.get(f"info/software", params={}, format=format)

    def info_species(self, format="json", division=None, hide_strain_info=None, strain_collection=None):
        """Lists all available species, their aliases, available adaptor groups and data release."""
        return self.get(f"info/species", params=dict(division=division, hide_strain_info=hide_strain_info, strain_collection=strain_collection), format=format)

    def info_variation(self, species, format="json", filter=None):
        """List the variation sources used in Ensembl for a species."""
        return self.get(f"info/variation/{species}", params=dict(filter=filter), format=format)

    def info_variation_consequence_types(self, format="json", rank=None):
        """Lists all variant consequence types."""
        return self.get(f"info/variation/consequence_types", params=dict(rank=rank), format=format)

    def info_variation_populations(self, species, population_name, format="json"):
        """List all individuals for a population from a species"""
        return self.get(f"info/variation/populations/{species}/{population_name}", params={}, format=format)

    def info_variation_species(self, species, format="json", filter=None):
        """List all populations for a species"""
        return self.get(f"info/variation/populations/{species}", params=dict(filter=filter), format=format)

    def ld(self, species, id, population_name, format="json", attribs=None, callback=None, d_prime=None, r2=None, window_size=None):
        """Computes and returns LD values between the given variant and all other variants in a window centered around the given variant.
        The window size is set to 500 kb."""
        return self.get(f"ld/{species}/{id}/{population_name}", params=dict(attribs=attribs, callback=callback, d_prime=d_prime, r2=r2, window_size=window_size), format=format)

    def ld_pairwise(self, species, id1, id2, format="json", d_prime=None, population_name=None, r2=None):
        """Computes and returns LD values between the given variants."""
        return self.get(f"ld/{species}/pairwise/{id1}/{id2}", params=dict(d_prime=d_prime, population_name=population_name, r2=r2), format=format)

    def ld_region(self, species, region, population_name, format="json", d_prime=None, r2=None):
        """Computes and returns LD values between all pairs of variants in the defined region."""
        return self.get(f"ld/{species}/region/{region}/{population_name}", params=dict(d_prime=d_prime, r2=r2), format=format)

    @singledispatchmethod
    def lookup_id(self, id: str, format="json", db_type=None, expand=None, mane=None, phenotypes=None, species=None, utr=None):
        """Find the species and database for a single identifier e.g. gene, transcript, protein"""
        return self.get(endpoint=f"lookup/id/{id}", params=dict(db_type=db_type, expand=expand, mane=mane, phenotypes=phenotypes,
                                                                species=species, utr=utr), format=format)

    @lookup_id.register
    def _(self, id: list, format="json", db_type=None, expand=None, mane=None, phenotypes=None, species=None, utr=None):
        """Find the species and database for several identifiers. IDs that are not found are returned with no data."""
        return self.post(endpoint=f"lookup/id", params=dict(db_type=db_type, expand=expand, mane=mane, phenotypes=phenotypes,
                                                            species=species, utr=utr), json={"ids": id}, format=format)

    @singledispatchmethod
    def lookup_symbol(self, symbol: str, format="json", species="human", expand=None):
        """Find the species and database for a symbol in a linked external database"""
        return self.get(f"lookup/symbol/{species}/{symbol}", params=dict(expand=expand), format=format)

    @lookup_symbol.register
    def _(self, symbol: list, species="human", format="json", expand=None):
        """Find the species and database for a set of symbols in a linked external database. Unknown symbols are omitted from the response."""
        return self.post(f"lookup/symbol/{species}", params=dict(expand=expand), json={"symbols": symbol}, format=format)

    def map_cdna(self, id, region, format="json", include_original_region=None, species=None):
        """Convert from cDNA coordinates to genomic coordinates. Output reflects forward orientation coordinates as returned from the Ensembl API."""
        return self.get(f"map/cdna/{id}/{region}", params=dict(include_original_region=include_original_region, species=species), format=format)

    def map_cds(self, id, region, format="json", include_original_region=None, species=None):
        """Convert from CDS coordinates to genomic coordinates. Output reflects forward orientation coordinates as returned from the Ensembl API."""
        return self.get(f"map/cds/{id}/{region}", params=dict(include_original_region=include_original_region, species=species), format=format)

    def map_assembly(self, species, asm_one, region, asm_two, format="json", coord_system=None, target_coord_system=None):
        """Convert the co-ordinates of one assembly to another"""
        return self.get(f"map/{species}/{asm_one}/{region}/{asm_two}", params=dict(coord_system=coord_system, target_coord_system=target_coord_system), format=format)

    def map_translation(self, id, region, format="json", species=None):
        """Convert from protein (translation) coordinates to genomic coordinates. Output reflects forward orientation coordinates as returned from the Ensembl API."""
        return self.get(f"map/translation/{id}/{region}", params=dict(species=species), format=format)

    def ontology_ancestors(self, id, format="json", ontology=None):
        """Reconstruct the entire ancestry of a term from is_a and part_of relationships"""
        return self.get(f"ontology/ancestors/{id}", params=dict(ontology=ontology), format=format)

    def ontology_ancestors_chart(self, id, format="json", ontology=None):
        """Reconstruct the entire ancestry of a term from is_a and part_of relationships"""
        return self.get(f"ontology/ancestors/chart/{id}", params=dict(ontology=ontology), format=format)

    def ontology_descendants(self, id, format="json", closest_term=None, ontology=None, subset=None, zero_distance=None):
        """Find all the terms descended from a given term. By default searches are conducted within the namespace of the given identifier"""
        return self.get(f"ontology/descendants/{id}", params=dict(closest_term=closest_term, ontology=ontology, subset=subset, zero_distance=zero_distance), format=format)

    def ontology_id(self, id, format="json", relation=None, simple=None):
        """Search for an ontological term by its namespaced identifier"""
        return self.get(f"ontology/id/{id}", params=dict(relation=relation, simple=simple), format=format)

    def ontology_name(self, name, format="json", ontology=None, relation=None, simple=None):
        """Search for a list of ontological terms by their name"""
        return self.get(f"ontology/name/{name}", params=dict(ontology=ontology, relation=relation, simple=simple), format=format)

    def taxonomy_classification(self, id, format="json"):
        """Return the taxonomic classification of a taxon node"""
        return self.get(f"taxonomy/classification/{id}", params={}, format=format)

    def taxonomy_id(self, id, format="json", simple=None):
        """Search for a taxonomic term by its identifier or name"""
        return self.get(f"taxonomy/id/{id}", params=dict(simple=simple), format=format)

    def taxonomy_name(self, name, format="json"):
        """Search for a taxonomic id by a non-scientific name"""
        return self.get(f"taxonomy/name/{name}", params={}, format=format)

    def overlap_id(self, id, format="json", feature=None, biotype=None, db_type=None, logic_name=None, misc_set=None, object_type=None, so_term=None, species=None,
                   species_set=None, variant_set=None):
        """Retrieves features (e.g. genes, transcripts, variants and more) that overlap a region defined by the given identifier."""
        return self.get(f"overlap/id/{id}", params=dict(feature=feature, biotype=biotype, db_type=db_type, logic_name=logic_name, misc_set=misc_set, object_type=object_type, so_term=so_term, species=species,
                                                        species_set=species_set, variant_set=variant_set), format=format)

    def overlap_region(self, species, region, format="json", feature=None, biotype=None, db_type=None, logic_name=None, misc_set=None, so_term=None, species_set=None,
                       trim_downstream=None, trim_upstream=None, variant_set=None):
        """Retrieves features (e.g. genes, transcripts, variants and more) that overlap a given region."""
        return self.get(f"overlap/region/{species}/{region}", params=dict(feature=feature, biotype=biotype, db_type=db_type, logic_name=logic_name, misc_set=misc_set, so_term=so_term,
                                                                          species_set=species_set, trim_downstream=trim_downstream, trim_upstream=trim_upstream, variant_set=variant_set), format=format)

    def overlap_translation(self, id, format="json", db_type=None, feature=None, so_term=None, species=None, type=None):
        """Retrieve features related to a specific Translation as described by its stable ID (e.g. domains, variants)."""
        return self.get(f"overlap/translation/{id}", params=dict(db_type=db_type, feature=feature, so_term=so_term, species=species, type=type), format=format)

    def phenotype_accession(self, species, accession, format="json", include_children=None, include_pubmed_id=None, include_review_status=None, source=None):
        """Return phenotype annotations for genomic features given a phenotype ontology accession"""
        return self.get(f"phenotype/accession/{species}/{accession}", params=dict(include_children=include_children, include_pubmed_id=include_pubmed_id, include_review_status=include_review_status, source=source), format=format)

    def phenotype_gene(self, species, gene, format="json", include_associated=None, include_overlap=None, include_pubmed_id=None, include_review_status=None, include_submitter=None, non_specified=None, trait=None, tumour=None):
        """Return phenotype annotations for a given gene."""
        return self.get(f"phenotype/gene/{species}/{gene}", params=dict(include_associated=include_associated, include_overlap=include_overlap, include_pubmed_id=include_pubmed_id, include_review_status=include_review_status, include_submitter=include_submitter, non_specified=non_specified,
                                                                        trait=trait, tumour=tumour), format=format)

    def phenotype_region(self, species, region, format="json", feature_type=None, include_pubmed_id=None, include_review_status=None, include_submitter=None, non_specified=None, only_phenotypes=None, trait=None, tumour=None):
        """Return phenotype annotations that overlap a given genomic region."""
        return self.get(f"phenotype/region/{species}/{region}", params=dict(feature_type=feature_type, include_pubmed_id=include_pubmed_id, include_review_status=include_review_status, include_submitter=include_submitter, non_specified=non_specified, only_phenotypes=only_phenotypes, trait=trait, tumour=tumour), format=format)

    def phenotype_term(self, species, term, format="json", include_children=None, include_pubmed_id=None, include_review_status=None, source=None):
        """Return phenotype annotations for genomic features given a phenotype ontology term"""
        return self.get(f"phenotype/term/{species}/{term}", params=dict(include_children=include_children, include_pubmed_id=include_pubmed_id, include_review_status=include_review_status, source=source), format=format)

    def regulatory_microarray_vendor(self, species, microarray, vendor, format="json"):
        """Returns information about a specific microarray"""
        return self.get(f"regulatory/species/{species}/microarray/{microarray}/vendor/{vendor}", params={}, format=format)

    def regulatory_species(self, species, format="json"):
        """Returns information about all epigenomes available for the given species"""
        return self.get(f"regulatory/species/{species}/epigenome", params={}, format=format)

    def species_binding_matrix(self, species, binding_matrix_stable_id, format="json", unit=None):
        """Return the specified binding matrix"""
        return self.get(f"species/{species}/binding_matrix/{binding_matrix_stable_id}", params=dict(unit=unit), format=format)

    def regulatory_microarray(self, species, format="json"):
        """Returns information about all microarrays available for the given species"""
        return self.get(f"regulatory/species/{species}/microarray", params={}, format=format)

    def regulatory_probe(self, species, microarray, probe, format="json", gene=None, transcripts=None):
        """Returns information about a specific probe from a microarray"""
        return self.get(f"regulatory/species/{species}/microarray/{microarray}/probe/{probe}", params=dict(gene=gene, transcripts=transcripts), format=format)

    def regulatory_probe_set(self, species, microarray, probe_set, format="json", gene=None, transcripts=None):
        """Returns information about a specific probe_set from a microarray"""
        return self.get(f"regulatory/species/{species}/microarray/{microarray}/probe/{probe_set}", params=dict(gene=gene, transcripts=transcripts), format=format)

    def regulatory_id(self, species, id, format="json", activity=None):
        """Returns a RegulatoryFeature given its stable ID (e.g. ENSR00000082023)"""
        return self.get(f"regulatory/species/{species}/id/{id}", params=dict(activity=activity), format=format)

    @singledispatchmethod
    def sequence_id(self, id: str, format="json", db_type=None, end=None, expand_3prime=None, expand_5prime=None, mask=None, mask_feature=None, multiple_sequences=None, object_type=None, species=None, start=None, type=None):
        """Request multiple types of sequence by stable identifier. Supports feature masking and expand options."""
        return self.get(f"sequence/id/{id}", params=dict(db_type=db_type, end=end, expand_3prime=expand_3prime, expand_5prime=expand_5prime, mask=mask, mask_feature=mask_feature, multiple_sequences=multiple_sequences, object_type=object_type, species=species, start=start, type=type), format=format)

    @sequence_id.register
    def _(self, id: list, format="json", db_type=None, end=None, expand_3prime=None, expand_5prime=None, mask=None, mask_feature=None, multiple_sequences=None, object_type=None, species=None, start=None, type=None):
        """Request multiple types of sequence by a stable identifier list."""
        return self.post(f"sequence/id", json={'ids': id}, params=dict(db_type=db_type, end=end, expand_3prime=expand_3prime, expand_5prime=expand_5prime, mask=mask, mask_feature=mask_feature, multiple_sequences=multiple_sequences, object_type=object_type, species=species, start=start, type=type), format=format)

    @singledispatchmethod
    def sequence_region(self, region: str, species, format="json", coord_system=None, coord_system_version=None, expand_3prime=None, expand_5prime=None, mask=None, mask_feature=None):
        """Returns the genomic sequence of the specified region of the given species. Supports feature masking and expand options."""
        return self.get(f"sequence/region/{species}/{region}", params=dict(coord_system=coord_system, coord_system_version=coord_system_version, expand_3prime=expand_3prime, expand_5prime=expand_5prime, mask=mask, mask_feature=mask_feature), format=format)

    @sequence_region.register
    def _(self, region: list, species, format="json", coord_system=None, coord_system_version=None, expand_3prime=None, expand_5prime=None, mask=None, mask_feature=None):
        """Request multiple types of sequence by a list of regions."""
        return self.post(f"sequence/region/{species}", json={'region': region}, params=dict(coord_system=coord_system, coord_system_version=coord_system_version, expand_3prime=expand_3prime, expand_5prime=expand_5prime, mask=mask, mask_feature=mask_feature), format=format)

    def transcript_haplotypes(self, species, id, format="json", aligned_sequences=None, samples=None, sequence=None):
        """Computes observed transcript haplotype sequences based on phased genotype data"""
        return self.get(f"transcript_haplotypes/{species}/{id}", params=dict(aligned_sequences=aligned_sequences, samples=samples, sequence=sequence), format=format)

    def ga4gh_beacon(self, format="json"):
        """Return Beacon information"""
        return self.get(f"ga4gh/beacon", params={}, format=format)

    @singledispatchmethod
    def ga4gh_beacon_query(self, query: str,  alternateBases=None, assemblyId=None, referenceBases=None, referenceName=None, start=None, end=None, endMax=None, endMin=None, startMax=None, startMin=None, variantType=None,
                           datasetIds=None, includeDatasetResponses=None, format="json"):
        """Return the Beacon response for allele information"""
        return self.get(f"ga4gh/beacon/{query}", params=dict(alternateBases=alternateBases, assemblyId=assemblyId, referenceBases=referenceBases, referenceName=referenceName, start=start, datasetIds=datasetIds, end=end, endMax=endMax, endMin=endMin, startMax=startMax, startMin=startMin, variantType=variantType, includeDatasetResponses=includeDatasetResponses), format=format)

    @ga4gh_beacon_query.register
    def _(self, query: list, format="json", alternateBases=None, assemblyId=None, referenceBases=None, referenceName=None, start=None, end=None, endMax=None, endMin=None, startMax=None, startMin=None, variantType=None,
          datasetIds=None, includeDatasetResponses=None):
        """Return the Beacon response for allele information"""
        return self.post(f"ga4gh/beacon", json={"query": query}, params=dict(alternateBases=alternateBases, assemblyId=assemblyId, referenceBases=referenceBases, referenceName=referenceName, start=start, datasetIds=datasetIds, end=end, endMax=endMax, endMin=endMin, startMax=startMax, startMin=startMin, variantType=variantType, includeDatasetResponses=includeDatasetResponses), format=format)

    def ga4gh_features(self, id, format="json"):
        """Return the GA4GH record for a specific sequence feature given its identifier"""
        return self.get(f"ga4gh/features/{id}", params={}, format=format)

    def ga4gh_features_search(self, end, referenceName, start, format="json", featureTypes=None, featureSetId=None, pageSize=None, pageToken=None, parentId=None):
        """Return a list of sequence annotation features in GA4GH format"""
        return self.post(f"ga4gh/features/search", json=dict(end=end, referenceName=referenceName, start=start, featureTypes=featureTypes, featureSetId=featureSetId, pageSize=pageSize, pageToken=pageToken, parentId=parentId), params={}, format=format)

    def ga4gh_callsets_search(self, variantSetId, format="json", **kwargs):
        """Return a list of sets of genotype calls for specific samples in GA4GH format"""
        return self.post(f"ga4gh/callsets/search", json=dict(variantSetId=variantSetId, **kwargs), params={}, format=format)

    def ga4gh_callsets(self, id, format="json"):
        """Return the GA4GH record for a specific CallSet given its identifier"""
        return self.get(f"ga4gh/callsets/{id}", params={}, format=format)

    def ga4gh_datasets_search(self, format="json", pageSize=None, pageToken=None):
        """Return a list of datasets in GA4GH format"""
        return self.post(f"ga4gh/datasets/search", json=dict(pageSize=pageSize, pageToken=pageToken), params={}, format=format)

    def ga4gh_datasets(self, id, format="json"):
        """Return the GA4GH record for a specific dataset given its identifier"""
        return self.get(f"ga4gh/datasets/{id}", params={}, format=format)

    def ga4gh_featuresets_search(self, datasetId, format="json", pageSize=None, pageToken=None):
        """Return a list of feature sets in GA4GH format"""
        return self.post(f"ga4gh/featuresets/search", json=dict(datasetId=datasetId, pageSize=pageSize, pageToken=pageToken), params={}, format=format)

    def ga4gh_featuresets(self, id, format="json"):
        """Return the GA4GH record for a specific featureSet given its identifier"""
        return self.get(f"ga4gh_featuresets/{id}", params={}, format=format)

    def ga4gh_variants(self, id, format="json"):
        """Return the GA4GH record for a specific variant given its identifier."""
        return self.get(f"ga4gh/variants/{id}", params={}, format=format)

    def ga4gh_variants_search(self, end, referenceName, start, variantSetId, format="json", callSetIds=None, pageSize=None, pageToken=None):
        """Return variant call information in GA4GH format for a region on a reference sequence"""
        return self.post(f"ga4gh/variants/search", json=dict(end=end, referenceName=referenceName, start=start, variantSetId=variantSetId, callSetIds=callSetIds, pageSize=pageSize, pageToken=pageToken), params={}, format=format)

    def ga4gh_variantannotations_search(self, variantAnnotationSetId, format="json", effects=None, end=None, pageSize=None, pageToken=None, referenceId=None, referenceName=None, start=None):
        """Return variant annotation information in GA4GH format for a region on a reference sequence"""
        return self.post(f"ga4gh/variantannotations/search", json=dict(variantAnnotationSetId=variantAnnotationSetId, effects=effects, end=end, pageSize=pageSize, pageToken=pageToken, referenceId=referenceId, referenceName=referenceName, start=start), params={}, format=format)

    def ga4gh_variantsets_search(self, datasetId, format="json", pageSize=None, pageToken=None):
        """Return a list of variant sets in GA4GH format"""
        return self.post(f"ga4gh/variantsets/search", json=dict(datasetId=datasetId, pageSize=pageSize, pageToken=pageToken), params={}, format=format)

    def ga4gh_variantsets(self, id, format="json"):
        """Return the GA4GH record for a specific VariantSet given its identifier"""
        return self.get(f"ga4gh/variantsets/{id}", params={}, format=format)

    def ga4gh_references_search(self, referenceSetId, format="json", accession=None, md5checksum=None, pageSize=None, pageToken=None):
        """Return a list of reference sequences in GA4GH format"""
        return self.post(f"ga4gh/references/search", json=dict(referenceSetId=referenceSetId, accession=accession, md5checksum=md5checksum, pageSize=pageSize, pageToken=pageToken), params={}, format=format)

    def ga4gh_references(self, id, format="json"):
        """Return data for a specific reference in GA4GH format by id"""
        return self.get(f"ga4gh/references/{id}", params={}, format=format)

    def ga4gh_referencesets_search(self, format="json", accession=None, pageSize=None, pageToken=None):
        """Return a list of reference sets in GA4GH format"""
        return self.post(f"ga4gh/referencesets/search", json=dict(accession=accession, pageSize=pageSize, pageToken=pageToken), params={}, format=format)

    def ga4gh_referencesets(self, id, format="json"):
        """Return data for a specific reference set in GA4GH format"""
        return self.get(f"ga4gh/referencesets/{id}", params={}, format=format)

    def ga4gh_variantannotationsets_search(self, variantSetId, format="json", pageSize=None, pageToken=None):
        """Return a list of annotation sets in GA4GH format"""
        return self.post(f"ga4gh/variantannotationsets/search", json=dict(variantSetId=variantSetId, pageSize=pageSize, pageToken=pageToken), params={}, format=format)

    def ga4gh_variantannotationsets(self, id, format="json"):
        """Return meta data for a specific annotation set in GA4GH format"""
        return self.get(f"ga4gh/variantannotationsets/{id}", params={}, format=format)

    @singledispatchmethod
    def variant_recoder(self, id: str, species='human', format="json", fields=None, var_synonyms=None, vcf_string=None):
        """Translate a variant identifier, HGVS notation or genomic SPDI notation to all possible variant IDs, HGVS and genomic SPDI"""
        return self.get(endpoint=f"variant_recoder/{species}/{id}", params=dict(fields=fields, var_synonyms=var_synonyms, vcf_string=vcf_string), format=format)

    @variant_recoder.register
    def _(self, id: list, species='human', format="json", fields=None, var_synonyms=None, vcf_string=None):
        """Translate a list of variant identifiers, HGVS notations or genomic SPDI notations to all possible variant IDs, HGVS and genomic SPDI"""
        return self.post(endpoint=f"variant_recoder/{species}", params=dict(fields=fields, var_synonyms=var_synonyms, vcf_string=vcf_string), json={"ids": id}, format=format)

    @singledispatchmethod
    def variation(self, id: str, species='human', format="json", pops=None, genotypes=None, genotyping_chips=None, phenotypes=None, population_genotypes=None):
        """Uses a variant identifier (e.g. rsID) to return the variation features including optional genotype, phenotype and population data"""
        return self.get(endpoint=f"variation/{species}/{id}", params=dict(pops=pops, genotypes=genotypes, genotyping_chips=genotyping_chips, phenotypes=phenotypes, population_genotypes=population_genotypes), format=format)

    @variation.register
    def _(self, id: list, species='human', format="json", pops=None, genotypes=None, genotyping_chips=None, phenotypes=None, population_genotypes=None):
        """Uses a list of variant identifiers (e.g. rsID) to return the variation features including optional genotype, phenotype and population data"""
        return self.post(endpoint=f"variation/{species}", params=dict(pops=pops, genotypes=genotypes, genotyping_chips=genotyping_chips, phenotypes=phenotypes, population_genotypes=population_genotypes), json={"ids": id}, format=format)

    def variation_pmcid(self, pmcid, species='human', format="json"):
        """Fetch variants by publication using PubMed Central reference number (PMCID)"""
        return self.get(endpoint=f"variation/{species}/pmcid/{pmcid}", params={}, format=format)

    def variation_pmid(self, pmid, species='human', format="json"):
        """Fetch variants by publication using PubMed reference number (PMID)"""
        return self.get(endpoint=f"variation/{species}/pmid/{pmid}", params={}, format=format)

    @singledispatchmethod
    def vep_hgvs(self, hgvs_notation: str, species='human', format="json", AncestralAllele=None, Blosum62=None, CADD=None, Conservation=None, DisGeNET=None, EVE=None, GO=None, GeneSplicer=None, IntAct=None, LoF=None, Mastermind=None, MaxEntScan=None,
                 NMD=None, Phenotypes=None, SpliceAI=None, UTRAnnotator=None, ambiguous_hgvs=None, appris=None, canonical=None, ccds=None, dbNSFP=None, dbscSNV=None, distance=None, domains=None, failed=None, hgvs=None, mane=None, merged=None, minimal=None,
                 mirna=None, mutfunc=None, numbers=None, protein=None, refseq=None, shift_3prime=None, shift_genomic=None, transcript_id=None, transcript_version=None, tsl=None, uniprot=None, variant_class=None, vcf_string=None, xref_refseq=None):
        """Fetch variant consequences based on a HGVS notation"""
        return self.get(endpoint=f"vep/{species}/hgvs/{hgvs_notation}", params=dict(AncestralAllele=AncestralAllele, Blosum62=Blosum62, CADD=CADD, Conservation=Conservation, DisGeNET=DisGeNET, EVE=EVE, GO=GO, GeneSplicer=GeneSplicer, IntAct=IntAct, LoF=LoF, Mastermind=Mastermind, MaxEntScan=MaxEntScan,
                                                                                    NMD=NMD, Phenotypes=Phenotypes, SpliceAI=SpliceAI, UTRAnnotator=UTRAnnotator, ambiguous_hgvs=ambiguous_hgvs, appris=appris, canonical=canonical, ccds=ccds, dbNSFP=dbNSFP, dbscSNV=dbscSNV, distance=distance, domains=domains, failed=failed, hgvs=hgvs, mane=mane, merged=merged, minimal=minimal,
                                                                                    mirna=mirna, mutfunc=mutfunc, numbers=numbers, protein=protein, refseq=refseq, shift_3prime=shift_3prime, shift_genomic=shift_genomic, transcript_id=transcript_id, transcript_version=transcript_version, tsl=tsl, uniprot=uniprot, variant_class=variant_class, vcf_string=vcf_string, xref_refseq=xref_refseq), format=format)

    @vep_hgvs.register
    def _(self, hgvs_notation: list, species='human', format="json", AncestralAllele=None, Blosum62=None, CADD=None, Conservation=None, DisGeNET=None, EVE=None, GO=None, GeneSplicer=None, IntAct=None, LoF=None, Mastermind=None, MaxEntScan=None,
          NMD=None, Phenotypes=None, SpliceAI=None, UTRAnnotator=None, ambiguous_hgvs=None, appris=None, canonical=None, ccds=None, dbNSFP=None, dbscSNV=None, distance=None, domains=None, failed=None, hgvs=None, mane=None, merged=None, minimal=None,
          mirna=None, mutfunc=None, numbers=None, protein=None, refseq=None, shift_3prime=None, shift_genomic=None, transcript_id=None, transcript_version=None, tsl=None, uniprot=None, variant_class=None, vcf_string=None, xref_refseq=None):
        """Fetch variant consequences for multiple HGVS notations"""
        return self.post(endpoint=f"vep/{species}/hgvs", params=dict(AncestralAllele=AncestralAllele, Blosum62=Blosum62, CADD=CADD, Conservation=Conservation, DisGeNET=DisGeNET, EVE=EVE, GO=GO, GeneSplicer=GeneSplicer, IntAct=IntAct, LoF=LoF, Mastermind=Mastermind, MaxEntScan=MaxEntScan,
                                                                     NMD=NMD, Phenotypes=Phenotypes, SpliceAI=SpliceAI, UTRAnnotator=UTRAnnotator, ambiguous_hgvs=ambiguous_hgvs, appris=appris, canonical=canonical, ccds=ccds, dbNSFP=dbNSFP, dbscSNV=dbscSNV, distance=distance, domains=domains, failed=failed, hgvs=hgvs, mane=mane, merged=merged, minimal=minimal,
                                                                     mirna=mirna, mutfunc=mutfunc, numbers=numbers, protein=protein, refseq=refseq, shift_3prime=shift_3prime, shift_genomic=shift_genomic, transcript_id=transcript_id, transcript_version=transcript_version, tsl=tsl, uniprot=uniprot, variant_class=variant_class, vcf_string=vcf_string, xref_refseq=xref_refseq), json={"hgvs_notations": hgvs_notation}, format=format)

    @singledispatchmethod
    def vep_id(self, id: str, species='human', format="json", AncestralAllele=None, Blosum62=None, CADD=None, Conservation=None, DisGeNET=None, EVE=None, GO=None, GeneSplicer=None, IntAct=None, LoF=None, Mastermind=None, MaxEntScan=None,
               NMD=None, Phenotypes=None, SpliceAI=None, UTRAnnotator=None, appris=None, canonical=None, ccds=None, dbNSFP=None, dbscSNV=None, distance=None, domains=None, failed=None, hgvs=None, mane=None, merged=None, minimal=None,
               mirna=None, mutfunc=None, numbers=None, protein=None, refseq=None, shift_3prime=None, shift_genomic=None, transcript_id=None, transcript_version=None, tsl=None, uniprot=None, variant_class=None, vcf_string=None, xref_refseq=None):
        """Fetch variant consequences based on a variant identifier"""
        return self.get(endpoint=f"vep/{species}/id/{id}", params=dict(AncestralAllele=AncestralAllele, Blosum62=Blosum62, CADD=CADD, Conservation=Conservation, DisGeNET=DisGeNET, EVE=EVE, GO=GO, GeneSplicer=GeneSplicer, IntAct=IntAct, LoF=LoF, Mastermind=Mastermind, MaxEntScan=MaxEntScan,
                                                                       NMD=NMD, Phenotypes=Phenotypes, SpliceAI=SpliceAI, UTRAnnotator=UTRAnnotator, appris=appris, canonical=canonical, ccds=ccds, dbNSFP=dbNSFP, dbscSNV=dbscSNV, distance=distance, domains=domains, failed=failed, hgvs=hgvs, mane=mane, merged=merged, minimal=minimal,
                                                                       mirna=mirna, mutfunc=mutfunc, numbers=numbers, protein=protein, refseq=refseq, shift_3prime=shift_3prime, shift_genomic=shift_genomic, transcript_id=transcript_id, transcript_version=transcript_version, tsl=tsl, uniprot=uniprot, variant_class=variant_class, vcf_string=vcf_string, xref_refseq=xref_refseq), format=format)

    @vep_id.register
    def _(self, id: list, species='human', format="json", AncestralAllele=None, Blosum62=None, CADD=None, Conservation=None, DisGeNET=None, EVE=None, GO=None, GeneSplicer=None, IntAct=None, LoF=None, Mastermind=None, MaxEntScan=None,
            NMD=None, Phenotypes=None, SpliceAI=None, UTRAnnotator=None, appris=None, canonical=None, ccds=None, dbNSFP=None, dbscSNV=None, distance=None, domains=None, failed=None, hgvs=None, mane=None, merged=None, minimal=None,
            mirna=None, mutfunc=None, numbers=None, protein=None, refseq=None, shift_3prime=None, shift_genomic=None, transcript_id=None, transcript_version=None, tsl=None, uniprot=None, variant_class=None, vcf_string=None, xref_refseq=None):
        """Fetch variant consequences for multiple ids"""
        return self.post(endpoint=f"vep/{species}/id", params=dict(AncestralAllele=AncestralAllele, Blosum62=Blosum62, CADD=CADD, Conservation=Conservation, DisGeNET=DisGeNET, EVE=EVE, GO=GO, GeneSplicer=GeneSplicer, IntAct=IntAct, LoF=LoF, Mastermind=Mastermind, MaxEntScan=MaxEntScan,
                                                                   NMD=NMD, Phenotypes=Phenotypes, SpliceAI=SpliceAI, UTRAnnotator=UTRAnnotator, appris=appris, canonical=canonical, ccds=ccds, dbNSFP=dbNSFP, dbscSNV=dbscSNV, distance=distance, domains=domains, failed=failed, hgvs=hgvs, mane=mane, merged=merged, minimal=minimal,
                                                                   mirna=mirna, mutfunc=mutfunc, numbers=numbers, protein=protein, refseq=refseq, shift_3prime=shift_3prime, shift_genomic=shift_genomic, transcript_id=transcript_id, transcript_version=transcript_version, tsl=tsl, uniprot=uniprot, variant_class=variant_class, vcf_string=vcf_string, xref_refseq=xref_refseq), json={"ids": id}, format=format)

    @singledispatchmethod
    def vep_region(self, region: str, allele, species='human', format="json", AncestralAllele=None, Blosum62=None, CADD=None, Conservation=None, DisGeNET=None, EVE=None, GO=None, GeneSplicer=None, IntAct=None, LoF=None, Mastermind=None, MaxEntScan=None,
                   NMD=None, Phenotypes=None, SpliceAI=None, UTRAnnotator=None, appris=None, canonical=None, ccds=None, dbNSFP=None, dbscSNV=None, distance=None, domains=None, failed=None, hgvs=None, mane=None, merged=None, minimal=None,
                   mirna=None, mutfunc=None, numbers=None, protein=None, refseq=None, shift_3prime=None, shift_genomic=None, transcript_id=None, transcript_version=None, tsl=None, uniprot=None, variant_class=None, vcf_string=None, xref_refseq=None):
        """Fetch variant consequences based on a region"""
        return self.get(endpoint=f"vep/{species}/region/{region}/{allele}", params=dict(AncestralAllele=AncestralAllele, Blosum62=Blosum62, CADD=CADD, Conservation=Conservation, DisGeNET=DisGeNET, EVE=EVE, GO=GO, GeneSplicer=GeneSplicer, IntAct=IntAct, LoF=LoF, Mastermind=Mastermind, MaxEntScan=MaxEntScan,
                                                                                        NMD=NMD, Phenotypes=Phenotypes, SpliceAI=SpliceAI, UTRAnnotator=UTRAnnotator, appris=appris, canonical=canonical, ccds=ccds, dbNSFP=dbNSFP, dbscSNV=dbscSNV, distance=distance, domains=domains, failed=failed, hgvs=hgvs, mane=mane, merged=merged, minimal=minimal,
                                                                                        mirna=mirna, mutfunc=mutfunc, numbers=numbers, protein=protein, refseq=refseq, shift_3prime=shift_3prime, shift_genomic=shift_genomic, transcript_id=transcript_id, transcript_version=transcript_version, tsl=tsl, uniprot=uniprot, variant_class=variant_class, vcf_string=vcf_string, xref_refseq=xref_refseq), format=format)

    @vep_region.register
    def _(self, region: list, species='human', format="json", AncestralAllele=None, Blosum62=None, CADD=None, Conservation=None, DisGeNET=None, EVE=None, GO=None, GeneSplicer=None, IntAct=None, LoF=None, Mastermind=None, MaxEntScan=None,
          NMD=None, Phenotypes=None, SpliceAI=None, UTRAnnotator=None, appris=None, canonical=None, ccds=None, dbNSFP=None, dbscSNV=None, distance=None, domains=None, failed=None, hgvs=None, mane=None, merged=None, minimal=None,
          mirna=None, mutfunc=None, numbers=None, protein=None, refseq=None, shift_3prime=None, shift_genomic=None, transcript_id=None, transcript_version=None, tsl=None, uniprot=None, variant_class=None, vcf_string=None, xref_refseq=None):
        """Fetch variant consequences for multiple regions"""
        return self.post(endpoint=f"vep/{species}/region", params=dict(AncestralAllele=AncestralAllele, Blosum62=Blosum62, CADD=CADD, Conservation=Conservation, DisGeNET=DisGeNET, EVE=EVE, GO=GO, GeneSplicer=GeneSplicer, IntAct=IntAct, LoF=LoF, Mastermind=Mastermind, MaxEntScan=MaxEntScan,
                                                                       NMD=NMD, Phenotypes=Phenotypes, SpliceAI=SpliceAI, UTRAnnotator=UTRAnnotator, appris=appris, canonical=canonical, ccds=ccds, dbNSFP=dbNSFP, dbscSNV=dbscSNV, distance=distance, domains=domains, failed=failed, hgvs=hgvs, mane=mane, merged=merged, minimal=minimal,
                                                                       mirna=mirna, mutfunc=mutfunc, numbers=numbers, protein=protein, refseq=refseq, shift_3prime=shift_3prime, shift_genomic=shift_genomic, transcript_id=transcript_id, transcript_version=transcript_version, tsl=tsl, uniprot=uniprot, variant_class=variant_class, vcf_string=vcf_string, xref_refseq=xref_refseq), json={"variants": region}, format=format)
