#include <sstream>
#include <algorithm>
#include <iostream>
#include <vector>

bool IsValidDNASequence(const std::string& sequence) {
	/* This function uses a range based for loop to make sure all characters are ATCG.
	No spaces or any other characters are allowed.*/
	for (auto x : sequence) {
		if ((x != 'A') && (x != 'T') && (x != 'C') && (x != 'G'))
			return false;
	}
	return true;

}

void GetReverseComplementSequence(const std::string& input, std::string* const output) {
	/* Using a stringstream, this function concenates the Reverse Complement of a DNA sequence
	in order to later find the RNA transcript.*/
	std::string reverse = input;
	std::reverse(reverse.begin(), reverse.end());
	std::istringstream iss(reverse);
	char x;
	while (iss >> x) {
		switch (x) {
		case('A'):
			*output += 'T';
			break;
		case('T'):
			*output += 'A';
			break;
		case('C'):
			*output += 'G';
			break;
		case('G'):
			*output += 'C';
			break;
		}
	}
}
std::string GetRNATranscript(const std::string& input) {
	/* Using the Reverse Complement (returned via a pointer), this function uses a stringstream to change
	any instance of T to a U for the (original) RNA transcript*/
	std::string reverse_complement = "";
	std::string* output = &reverse_complement;
	std::string rna_transcript = "";
	char x;
	GetReverseComplementSequence(input, output);
	std::istringstream iss(*output);
	while (iss >> x) {
		if (x == 'T')
			rna_transcript += 'U';
		else
			rna_transcript += x;
	}
	return rna_transcript;
}



std::vector<std::vector<std::string>> GetReadingFramesAsCodons(const std::string& input) {
	/* This function finds the original and antiparallel RNA transcripts of a DNA transcript and
	uses the RNA transcripts to find codons - string triplets. Those are then put into a 2D vector.*/
	std::string offset_transcript, codon, old_offset;
	std::string offset = input;
	std::vector<std::vector<std::string>> main_vector;
	std::vector<std::string> inner_vector;
	char character;
	for (int i = 0; i <= 2; i++) {  // For loop runs 3 times for the 3 offsets.
		offset_transcript = GetRNATranscript(offset);
		if (offset_transcript != old_offset)
			codon = ""; // Clear the codon string to add a new one
		std::istringstream iss(offset_transcript);
		while (iss >> character) {
			codon += character;
			if (codon.length() == 3) {
				inner_vector.push_back(codon);
				codon = "";
				offset_transcript.erase(0, 2); // Removes the first element of the transcript - this creates the offsets.
			}
		}
		main_vector.push_back(inner_vector);
		old_offset = offset_transcript;
		offset.pop_back();
		inner_vector.clear(); // Empties the vector to allow others to be pushed back
	}
	std::string antiparallel = input;
	std::replace(antiparallel.begin(), antiparallel.end(), 'T', 'U'); // Quickly replaces Ts with Us for the antiparallel strand
	offset = antiparallel;
	for (int i = 0; i <= 2; i++) { // Same process as above - I'm sorry! I struggled with this part.
		if (offset != old_offset)
			codon = "";
		std::istringstream iss(offset);
		while (iss >> character) {
			codon += character;
			if (codon.length() == 3) {
				inner_vector.push_back(codon);
				codon = "";
			}
		}
		main_vector.push_back(inner_vector);
		old_offset = offset;
		offset.erase(0, 1);
		inner_vector.clear();
	}

	return main_vector;
}


std::string Translate(const std::vector<std::string>& codon_sequence) {
	std::vector<std::string> codon_list = { "GCU", "GCC", "GCA", "GCG", "CGU", "CGC", "CGA", "CGG", "AGA", "AGG",
		"AAU", "AAC", "GAU", "GAC", "UGU", "UGC", "CAA", "CAG", "GAA", "GAG",
		"GGU", "GGC", "GGA", "GGG", "CAU", "CAC", "AUU", "AUC", "AUA", "UUA",
		"UUG", "CUU", "CUC", "CUA", "CUG", "AAA", "AAG", "AUG", "UUU", "UUC",
		"CCU", "CCC", "CCA", "CCG", "UCU", "UCC", "UCA", "UCG", "AGU", "AGC",
		"ACU", "ACC", "ACA", "ACG", "UGG", "UAU", "UAC", "GUU", "GUC", "GUA",
		"GUG", "UAG", "UGA", "UAA" };
	std::vector<std::string> corresponding_codons = { "A", "A", "A", "A", "R", "R", "R", "R", "R", "R", "N", "N", "D", "D",
		"C", "C", "Q", "Q", "E", "E", "G", "G", "G", "G", "H", "H", "I", "I",
		"I", "L", "L", "L", "L", "L", "L", "K", "K", "M", "F", "F", "P", "P",
		"P", "P", "S", "S", "S", "S", "S", "S", "T", "T", "T", "T", "W", "Y",
		"Y", "V", "V", "V", "V", "*", "*", "*" };
	std::string translated_codon;

	for (auto codon : codon_sequence) {
		/* This for loop uses iterators in order to determine where a certain codon is located in the above vector.
		The same index is then used to find the corresponding codon in the (fittingly named) corresponding_codons vector
		Credit to Nitash on Piazza - see post https://piazza.com/class/k4wvjqt9b2x72k?cid=786. Thanks!*/
		auto codon_location = std::find(std::begin(codon_list), std::end(codon_list), codon);
		auto distance = std::distance(std::begin(codon_list), codon_location);
		translated_codon += corresponding_codons[distance];
	}
	return translated_codon;
}

std::string GetLongestOpenReadingFrame(const std::string& DNA_sequence) {
	std::vector<std::vector<std::string>> codon_vector = GetReadingFramesAsCodons(DNA_sequence);
	bool valid_reading_frame = false;
	int highest_count = 0;
	int count = 0;
	std::string longest_frame, current_frame;
	for (auto inner_vector : codon_vector) { 
		std::string translated_sequence = Translate(inner_vector); 
		for (auto element : translated_sequence) {
			if (element == 'M' && count == 0) { // Start of reading frame
				valid_reading_frame = true;
				current_frame += element; // Add the current character to the frame that is being read
			}
			else if (element != '*' && valid_reading_frame == true) { // Continue to count in a reading frame
				current_frame += element; 
				count++;
			}
			else if (element == '*' && valid_reading_frame == true) { // End reading frame
				current_frame += element;
				if (count > highest_count) { 
					highest_count = count;
					longest_frame = current_frame; // Copy the current frame if it is the longest, discard it otherwise
				}
				valid_reading_frame = false; // "Re-initialize" variables after stop codon is encountered
				count = 0;
				current_frame = "";

			}
		}
		valid_reading_frame = false; // "Re-initialize variables at end of sequence (or line)
		count = 0;
		current_frame = "";
	}
	return longest_frame;
}