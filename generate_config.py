import re

def generate_config():
    """
    Reads constants from assembler/constants.py and generates a config.ini file.
    """
    config_content = "[DEFAULT]\n"
    with open('assembler/constants.py', 'r') as f:
        for line in f:
            line = line.strip()
            # Ignore empty lines and comments
            if not line or line.startswith('#'):
                continue
            
            # Use regex to find variable assignments
            match = re.match(r"^(\w+)\s*=\s*(.*)", line)
            if match:
                key = match.group(1)
                value = match.group(2).strip()
                
                # Remove any inline comments from the value
                value = value.split('#')[0].strip()
                
                config_content += f"{key} = {value}\n"

    with open('config.ini', 'w') as f:
        f.write(config_content)
    
    print("config.ini has been generated successfully.")

if __name__ == "__main__":
    generate_config()
