#!/usr/bin/env python3

import argparse
import logging
import json
import os
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import re
import sys
import os
from datetime import datetime
import multiprocessing

import iso8601

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

ACTUALLY_MISSING = set([s.lower() for s in [
    'missing','not applicable','NA','Missing','Not collected','not provided','Missing: Not provided','', 'uncalculated','not applicable','no applicable','unspecified','restricted access']])


def validate_lat_lon(lat, lon):
    if lat >= -90 and lat <= 90 and lon >= -180 and lon <= 180:
        return True
    return False


def parse_json_to_lat_lon_dict(j):
    lat_long_dict = {}
	
    biosample_keys_for_lat_long_parsing = [
        'lat_lon_sam',
        'geographic_location__latitude__sam',
        'geographic_location__longitude__sam',
        'latitude_start_sam',
        'longitude_start_sam',
        'latitude_sam',
        'longitude_sam',
        'sampling_event__latitude__start_sam',
        'sampling_event__longitude__start_sam',
        'geographic_location__latitude_and_longitude__sam',
        'lat_lon_sam_s_dpl34'
    ]
	
    for attr in j['attributes']:
        if attr['k'] in biosample_keys_for_lat_long_parsing and attr['v'].lower() not in ACTUALLY_MISSING:
            lat_long_dict[attr['k']] = attr['v']
			
    return lat_long_dict


def parse_lat_lon_sam(lat_lon_sam):
    lat_lon_regex = re.compile(r'^([0-9.-]+) {0,1}([NSns]),{0,1} {0,1}([0-9.-]+) ([EWew])$')
    matches = lat_lon_regex.match(lat_lon_sam)
    if matches is None:
        logging.warning("Unexpected lat_lon_sam value: %s" % lat_lon_sam)
        return False, False
    try:
        lat = float(matches.group(1))
    except ValueError:
        logging.warning("Unexpected lat_lon_sam value: %s" % lat_lon_sam)
        return False, False
    if matches.group(2) in ['S','s']:
        lat = -lat
    try:
        lon = float(matches.group(3))
    except ValueError:
        logging.warning("Unexpected lat_lon_sam value: %s" % lat_lon_sam)
        return False, False
    if matches.group(4) in ['W','w']:
        lon = -lon
    if validate_lat_lon(lat, lon):
        return True, (lat, lon)
    else:
        logging.warning("Unvalidated lat_lon_sam value: %s" % lat_lon_sam)
        return False, False


def degrees_minutes_to_decimal(degrees, minutes, seconds=0):
    if seconds == '':
        seconds = 0
    return float(degrees) + float(minutes) / 60.0 + float(seconds) / 3600.0


def parse_two_part_lat_lon(sample_name, lat_input, lon_input):
    geoloc_lat_regex = re.compile(r'^([0-9.-]+)[°\?]{0,1} {0,1}([NSns])$')
    geoloc_lon_regex = re.compile(r'^([0-9.-]+)[°\?]{0,1} {0,1}([EWew])$')
    # e.g. S 12°37.707′
    sexigesimal_regex_lat1 = re.compile(r'^([NSns]) ([0-9]+)[°\?]([0-9.]+)[\'′]$')
    sexigesimal_regex_lon1 = re.compile(r'^([EWew]) ([0-9]+)[°\?]([0-9.]+)[\'′]$')
    # e.g. 52?09'50.8N
    sexigesimal_regex_lat2 = re.compile(r'^([0-9]+)[°\?]([0-9.]+)[\'′]([0-9.]*)([NSns])$')
    sexigesimal_regex_lon2 = re.compile(r'^([0-9]+)[°\?]([0-9.]+)[\'′]([0-9.]*)([EWew])$')
    # e.g.  ERR2824916: ["S33°28'21.68", "O70°38'50.06"] -> Actually that one is fail
    sexigesimal_regex_lat3 = re.compile(r'^([NSns])([0-9]+)[°\?]([0-9.]+)[\'′]([0-9.]*)$')
    sexigesimal_regex_lon3 = re.compile(r'^([EWew])([0-9]+)[°\?]([0-9.]+)[\'′]([0-9.]*)$')
    
    # e.g. ERR5866854: ["01° 30.464'N", "075° 55.202'W"] or ERR1956802: ["44°26'4272 N", "11°28'3540 E"] or ERR2044635: ["44°26'4272 N", "11°28'3540 E"]
    sexigesimal_regex_lat4 = re.compile(r'^([0-9]+)[°\?] ?([0-9]+(?:[.\'][0-9]+)?)[\'′ ]*([NSns])$')
    sexigesimal_regex_lon4 = re.compile(r'^([0-9]+)[°\?] ?([0-9]+(?:[.\'][0-9]+)?)[\'′ ]*([EWew])$')
    
    # e.g. ERR2271236: ["41° 53' 30 N", "12° 30' 40 E"]
    sexigesimal_regex_lat5 = re.compile(r'^([0-9]+)[°\?] ?([0-9]+)[\'′] ?([0-9]+) ?([NSns])$')
    sexigesimal_regex_lon5 = re.compile(r'^([0-9]+)[°\?] ?([0-9]+)[\'′] ?([0-9]+) ?([EWew])$')
    
    # e.g. ERR5173566: ['N 43.8047886', 'E 15.9637432']
    geoloc_lat_regex2 = re.compile(r'^([NSns]) {0,1}([0-9.-]+)[°\?]{0,1}$')
    geoloc_lon_regex2 = re.compile(r'^([EWew]) {0,1}([0-9.-]+)[°\?]{0,1}$')
    
        # Function to determine if a string contains latitude or longitude indicator
    def is_latitude(part):
        return bool(re.search(r'[NSns]', part))
    
    def is_longitude(part):
        return bool(re.search(r'[EWew]', part))
    
    # Check if lat_input and lon_input need to be swapped
    if is_latitude(lon_input) and is_longitude(lat_input):
        lat_input, lon_input = lon_input, lat_input
        
	
    try:
        lat = float(lat_input)
        lon = float(lon_input)
    except ValueError:
        # Try comma to dot and then convert
        try:
            lat = float(lat_input.replace(',','.'))
            lon = float(lon_input.replace(',','.'))
        except ValueError:
            matches_lat = geoloc_lat_regex.match(lat_input)
            matches_lon = geoloc_lon_regex.match(lon_input)
            if matches_lat is not None and matches_lon is not None:
                try:
                    lat = float(matches_lat.group(1))
                    lon = float(matches_lon.group(1))
                    if matches_lat.group(2) in ['S','s']:
                        lat = -lat
                    if matches_lon.group(2) in ['W','w']:
                        lon = -lon
                except ValueError:
                    logging.warning("Unexpected (type 1) 2 part value for %s: %s" % (sample_name, [lat_input, lon_input]))
                    return False, False
            else:
                # Try sexigesimal
                matches_lat1 = sexigesimal_regex_lat1.match(lat_input)
                matches_lon1 = sexigesimal_regex_lon1.match(lon_input)
                if matches_lat1 is not None and matches_lon1 is not None:
                    try:
                        lat = degrees_minutes_to_decimal(matches_lat1.group(2), matches_lat1.group(3))
                        lon = degrees_minutes_to_decimal(matches_lon1.group(2), matches_lon1.group(3))
						
                        if matches_lat1.group(1) in ['S','s']: lat = -lat
                        if matches_lon1.group(1) in ['W','w']: lon = -lon
                    except ValueError:
                        logging.warning("Unexpected (type 2) 2 part value for %s: %s" % (sample_name, [lat_input, lon_input]))
                        return False, False
                else:
                    matches_lat2 = sexigesimal_regex_lat2.match(lat_input)
                    matches_lon2 = sexigesimal_regex_lon2.match(lon_input)
                    if matches_lat2 is not None and matches_lon2 is not None:
                        try:
                            lat = degrees_minutes_to_decimal(matches_lat2.group(1), matches_lat2.group(2), matches_lat2.group(3))
                            lon = degrees_minutes_to_decimal(matches_lon2.group(1), matches_lon2.group(2), matches_lon2.group(3))
							
                            if matches_lat2.group(4) in ['S','s']: lat = -lat
                            if matches_lon2.group(4) in ['W','w']: lon = -lon
                        except ValueError:
                            logging.warning("Unexpected (type 4) 2 part value for %s: %s" % (sample_name, [lat_input, lon_input]))
                            return False, False
                    else:
                        matches_lat3 = sexigesimal_regex_lat3.match(lat_input)
                        matches_lon3 = sexigesimal_regex_lon3.match(lon_input)
                        if matches_lat3 is not None and matches_lon3 is not None:
                            try:
                                lat = degrees_minutes_to_decimal(matches_lat3.group(2), matches_lat3.group(3), matches_lat3.group(4))
                                lon = degrees_minutes_to_decimal(matches_lon3.group(2), matches_lon3.group(3), matches_lon3.group(4))
								
                                if matches_lat3.group(1) in ['S','s']: lat = -lat
                                if matches_lon3.group(1) in ['W','w']: lon = -lon
                            except ValueError:
                                logging.warning("Unexpected (type 5) 2 part value for %s: %s" % (sample_name, [lat_input, lon_input]))
                                return False, False
                        else:
                            matches_lat4 = geoloc_lat_regex2.match(lat_input)
                            matches_lon4 = geoloc_lon_regex2.match(lon_input)
                            if matches_lat4 is not None and matches_lon4 is not None:
                                try:
                                    lat = float(matches_lat4.group(2))
                                    lon = float(matches_lon4.group(2))
                                    if matches_lat4.group(1) in ['S','s']: lat = -lat
                                    if matches_lon4.group(1) in ['W','w']: lon = -lon
                                except ValueError:
                                    logging.warning("Unexpected (type 6) 2 part value for %s: %s" % (sample_name, [lat_input, lon_input]))
                                    return False, False
                            else:
                                matches_lat5 = sexigesimal_regex_lat4.match(lat_input)
                                matches_lon5 = sexigesimal_regex_lon4.match(lon_input)
                                if matches_lat5 is not None and matches_lon5 is not None:
                                    try:
                                        lat_parts = re.split('[°\'\\.]', matches_lat5.group(2).replace("'", ".").replace(" ", ""))
                                        lon_parts = re.split('[°\'\\.]', matches_lon5.group(2).replace("'", ".").replace(" ", ""))
                                        lat = degrees_minutes_to_decimal(matches_lat5.group(1), lat_parts[0], lat_parts[1] if len(lat_parts) > 1 else 0)
                                        lon = degrees_minutes_to_decimal(matches_lon5.group(1), lon_parts[0], lon_parts[1] if len(lon_parts) > 1 else 0)
                                        if matches_lat5.group(3) in ['S', 's']: lat = -lat
                                        if matches_lon5.group(3) in ['W', 'w']: lon = -lon
                                    except ValueError:
                                        logging.warning("Unexpected (type 6) 2 part value for %s: %s" % (sample_name, [lat_input, lon_input]))
                                        return False, False
                                else:
                                    matches_lat6 = sexigesimal_regex_lat5.match(lat_input)
                                    matches_lon6 = sexigesimal_regex_lon5.match(lon_input)
                                    if matches_lat6 is not None and matches_lon6 is not None:
                                        try:
                                            lat = degrees_minutes_to_decimal(matches_lat6.group(1), matches_lat6.group(2), matches_lat6.group(3))
                                            lon = degrees_minutes_to_decimal(matches_lon6.group(1), matches_lon6.group(2), matches_lon6.group(3))
                                            if matches_lat6.group(4) in ['S', 's']: lat = -lat
                                            if matches_lon6.group(4) in ['W', 'w']: lon = -lon
                                        except ValueError:
                                            logging.warning("Unexpected (type 7) 2 part value for %s: %s" % (sample_name, [lat_input, lon_input]))
                                            return False, False
                                    else:
                                        logging.warning("Unexpected (no regex match) 2 part value for %s: %s" % (sample_name, [lat_input, lon_input]))
                                        return False, False
								
    if lat and lon and validate_lat_lon(lat, lon):
        return True, (lat, lon)
    else:
        logging.warning("Unvalidated 2-part value: %s / %s" % (lat_input, lon_input))
        return False, False


def parse_temp(acc, temp):
    if temp.lower() in ACTUALLY_MISSING:
        return ''
		
    endings = ['°C','C','°c','c',' celcius']
    temp2 = temp
    for e in endings:
        if temp2.endswith(e):
            temp2 = temp2[:-len(e)]
			
    try:
        f = float(temp2)
        if f > 100: # SRR12345782 has 255C
            logging.warning("Too big temperature value for %s: %s from %s" % (acc, temp2, temp))
            return ''
        return f
    except ValueError:
        try:
            return float(temp2.replace(',','.'))
        except ValueError:
            logging.warning("Unexpected temperature value for %s: %s from %s" % (acc, temp2, temp))
            return ''


def parse_depth(acc, depth):
    if depth.lower() in ACTUALLY_MISSING:
        return ''
		
    depth2 = depth
	
    cm_endings = ['cm','centimeters','centimeter','centimetres','centimetre']
    is_cm = False # Assume metres
    for e in cm_endings:
        if depth2.endswith(e):
            depth2 = depth2[:-len(e)]
            is_cm = True
			
    if not is_cm:
        endings = ['m','meters','meter','metres','metre']
        for e in endings:
            if depth2.endswith(e):
                depth2 = depth2[:-len(e)]
				
    try:
        if is_cm:
            return float(depth2) / 100.0
        else:
            return float(depth2)
    except ValueError:
        logging.warning("Unexpected depth value for %s: %s from %s" % (acc, depth2, depth))
        return ''


def parse_date(acc, date):
    if date.lower() in ACTUALLY_MISSING:
        return ''
		
    try:
        d = iso8601.parse_date(date)
    except iso8601.ParseError:
        logging.warning("Unexpected date value for %s: %s" % (acc, date))
        return ['']*2
		
    if d.year < 1990 or d.year > datetime.now().year:
        logging.warning("Unexpected year value for %s: %s" % (acc, date))
    else:
        return [d.year, d.month]
		
    logging.warning("Unexpected date value for %s: %s" % (acc, date))
    return ['']*2



def parse_lat_lons(acc, lat_long_dict):
    single_part_keys = [
        'lat_lon_sam',
        'geographic_location__latitude_and_longitude__sam',
        'lat_lon_sam_s_dpl34'
    ]
    two_part_keys = [
        ['geographic_location__latitude__sam',
        'geographic_location__longitude__sam'],
        ['latitude_start_sam',
        'longitude_start_sam'],
        ['latitude_sam',
        'longitude_sam'],
        ['sampling_event__latitude__start_sam',
        'sampling_event__longitude__start_sam'],
    ]
	
    to_return = []
	
    got_a_lat_lon = False
    for k in single_part_keys:
        if k in lat_long_dict:
            returned = parse_lat_lon_sam(lat_long_dict[k])
            if returned is not None: # Happens when validation fails
                got_a_lat_lon, lat_lon = returned
                if got_a_lat_lon:
                    to_return.append(lat_lon[0])
                    to_return.append(lat_lon[1])
                    break
					
    if not got_a_lat_lon:
        for k in two_part_keys:
            if k[0] in lat_long_dict and k[1] in lat_long_dict:
                got_a_lat_lon, lat_lon = parse_two_part_lat_lon(acc, lat_long_dict[k[0]], lat_long_dict[k[1]])
                if got_a_lat_lon:
                    to_return.append(lat_lon[0])
                    to_return.append(lat_lon[1])
                    break
					
    if not got_a_lat_lon:
        to_return.append(None)
        to_return.append(None)
		
    return to_return

# Special parsing methods
parsing_hash = {
    'depth_sam': parse_depth,
    'temperature_sam': parse_temp,
    'collection_date_sam': parse_date,
}
# For when the number of fields returned by parsing != 1
non_standard_output_field_names_hash = {
    'collection_date_sam': ['collection_year', 'collection_month'],
}




# Function to process each JSON file and return rows and keys
def process_json_file(json_file):
    rows = []  # List to hold each row as a dictionary
    
    # Parse each line in the JSON file
    with open(json_file, 'r') as f:
        for line in f:
            j = json.loads(line)
            
            # Create a dictionary for this row
            row = j.copy()
            
            row.pop('attributes', None)
            
            # Handle lat/long if necessary
            lat_lon = parse_lat_lons(j.get('acc', ''), parse_json_to_lat_lon_dict(j))
            row['latitude'], row['longitude'] = lat_lon if lat_lon != [None, None] else ('', '')
            
            # Create the attributes_dict for all keys in 'attributes'
            attributes_dict = {attr['k']: attr['v'] for attr in j.get('attributes', []) if 'k' in attr and 'v' in attr}
            
            # Now update the attributes_dict for special parsing cases
            for key, value in attributes_dict.items():
                if value.lower() not in ACTUALLY_MISSING:
                    if key in parsing_hash:
                        parsed_value = parsing_hash[key](j['acc'], value)
                        if key in non_standard_output_field_names_hash:
                            # Handle non-standard fields that return multiple values
                            field_names = non_standard_output_field_names_hash[key]
                            for field_name, field_value in zip(field_names, parsed_value):
                                row[field_name] = field_value
                        else:
                            row[key] = parsed_value
                    else:
                        # No special parsing, just keep the original value
                        row[key] = value
            
            rows.append(row)
    
    return rows

# Define the columns you want to reorder
reordered_columns = ['acc','bioproject','experiment','sample_name',
                     'sample_acc','biosample','assay_type','geo_loc_name_country_calc',
                     'geo_loc_name_country_continent_calc','geo_loc_name_sam','organism',
                     'latitude','longitude','bases','bytes']


# Main function to run the script
def main(json_dir, output_file, threads, debug, quiet):
    # Setup logging
    loglevel = logging.DEBUG if debug else logging.ERROR if quiet else logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    
    # Collect all JSON files in the directory
    json_files = [os.path.join(json_dir, f) for f in os.listdir(json_dir) if f.endswith('.json')]
    
    # Make sure there are JSON files in the directory
    if not json_files:
        logging.error("No JSON files found in the directory.")
        return
    
    # Use multiprocessing to process files in parallel
    with multiprocessing.Pool(processes=threads) as pool:
        results = pool.map(process_json_file, json_files)
       
        # Combine all rows from all files
    all_rows = [row for result in results for row in result]
    
    # Directly convert the list of dictionaries to a DataFrame (pandas handles missing keys)
    combined_df = pd.DataFrame(all_rows)
    
    # Create a new column order
    new_order = reordered_columns + [col for col in combined_df.columns if col not in reordered_columns]
    combined_df = combined_df[new_order]
    
    # Save the DataFrame to CSV
    
    combined_df.to_csv(output_file, index=False, sep="\t")
    logging.info(f"Data successfully saved to {output_file}")
   

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process multiple JSON files and combine into a TSV.')
    parser.add_argument('--json-dir', required=True, help='Directory containing JSON files.')
    parser.add_argument('--output-file', required=True, help='Output TSV file.')
    parser.add_argument('--threads', type=int, default=4, help='Number of parallel threads to use.')
    parser.add_argument('--debug', action='store_true', help='Enable debug output.')
    parser.add_argument('--quiet', action='store_true', help='Suppress non-error output.')
    
    args = parser.parse_args()
    
    main(json_dir=args.json_dir, output_file=args.output_file, threads=args.threads, debug=args.debug, quiet=args.quiet)



