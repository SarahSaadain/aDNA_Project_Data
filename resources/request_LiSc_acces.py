import requests
from bs4 import BeautifulSoup
import argparse
 
# This script demonstrates how to perform a POST request with a CSRF token.
# It is a general example and is not intended for any specific website.
# The URL and form fields are placeholders.
 
# Function to get the CSRF token from a specified URL
def get_csrf_token(session, url):
    """
    Fetches the CSRF token from a web page.
 
    Args:
        session (requests.Session): The session object to use for the request.
        url (str): The URL of the page containing the form.
 
    Returns:
        tuple: A tuple containing the CSRF token and the session cookie.
    """
    try:
        # Perform a GET request to the URL to get the page content and cookies
        response = session.get(url)
        response.raise_for_status() # Raise an exception for bad status codes
 
        # Parse the HTML content using BeautifulSoup
        soup = BeautifulSoup(response.text, 'html.parser')
 
        # Find the input field with the name 'csrfmiddlewaretoken'
        csrf_token_input = soup.find('input', {'name': 'csrfmiddlewaretoken'})
 
        if csrf_token_input:
            csrf_token = csrf_token_input.get('value')
            print(f"✅ CSRF token found: {csrf_token}")
            return csrf_token, session.cookies
        else:
            print("❌ CSRF token not found on the page.")
            return None, None
 
    except requests.exceptions.RequestException as e:
        print(f"❌ Error during GET request: {e}")
        return None, None
 
# Function to perform the POST request
def perform_post_request(session, url, payload_token, cookies, data):
    """
    Performs a POST request with the given data and CSRF token.
 
    Args:
        session (requests.Session): The session object to use.
        url (str): The URL to send the POST request to.
        token (str): The CSRF token to include in the form data.
        cookies (requests.cookies.RequestsCookieJar): The cookies from the GET request.
        data (dict): The data payload to send.
    """
    # Add the CSRF token to the data payload
    data['csrfmiddlewaretoken'] = payload_token
 
    #print("Payload:", data)
    #print("Cookies:", cookies)
 
    # Define the headers for the POST request
    headers = {
        'Referer': 'https://lisc.univie.ac.at/firewall/',
        'Content-Type': 'application/x-www-form-urlencoded'
    }
 
    try:
        print(f"🔍 Sending POST request to {url} with data: {data}")
 
        # Perform the POST request
        post_response = session.post(url, data=data, headers=headers, cookies=cookies)
        post_response.raise_for_status() # Raise an exception for bad status codes
 
        print(f"✅ POST request successful!")
        print(f"➡️  Status Code: {post_response.status_code}")
        #print("➡️  Response Content:")
        #print(post_response.text[:200] + '...' if len(post_response.text) > 200 else post_response.text)
 
        print(("✅ A mail should be sent to the user shortly. Please check your inbox and confirm the activation."))
 
    except requests.exceptions.RequestException as e:
        print(f"❌ Error during POST request: {e}")
        print(f"➡️  Status Code: {e.response.status_code}" if e.response else "N/A")
        print("➡️  Response Content:")
        print(e.response.text if e.response else "N/A")
 
# Function to get the public IP address
def get_public_ip():
    """
    Fetches the current public IP address from a public API.
 
    Returns:
        str: The public IP address as a string, or None if an error occurs.
    """
    try:
        # Use a public service that returns the IP address in plain text
        response = requests.get('https://api.ipify.org')
        response.raise_for_status() # Raise an exception for bad status codes
        print(f"✅ Public IP address retrieved: {response.text}")
        return response.text
    except requests.exceptions.RequestException as e:
        print(f"❌ Error getting public IP: {e}")
        return None
 
def main():
    """
    Main function to run the script.
    """
    # Setup command line arguments
    parser = argparse.ArgumentParser(description="Demonstrates a simple POST request with CSRF token handling.")
    parser.add_argument('--username', required=True, help='The username to use for the request.')
    parser.add_argument('--ip', help='The IP address to use for the request.')
    args = parser.parse_args()
 
    # The target URL for the GET and POST requests
    target_url = "https://lisc.univie.ac.at/firewall/" # Placeholder URL, change this to your target
    post_url = "https://lisc.univie.ac.at/firewall/request_activate"
 
    # Start a new session to persist cookies
    with requests.Session() as s:
        print(f"🔍 Getting CSRF token from {target_url}...")
        csrf_token, cookies = get_csrf_token(s, target_url)
 
        if not csrf_token:
            print("❌ Unable to proceed without a CSRF token.")
            return
 
        if args.ip:
            print(f"✅ Using provided IP address: {args.ip}")
        else:
            print("🔍 No IP address provided, fetching public IP...")
            args.ip = get_public_ip()
            if not args.ip:
                print("❌ Unable to fetch public IP, exiting.")
                return
 
        # Replace with your specific data payload
        # This is the data that will be sent in the POST request body
        post_data = {
            'ip': args.ip,
            'username': args.username,
            'totp': ''  # Replace with your TOTP logic
        }
 
        print(f"\n🚀 Requesting IP access for {args.ip} and user {args.username}...")
        perform_post_request(s, post_url, csrf_token, cookies, post_data)
 
if __name__ == "__main__":
    main()
 
