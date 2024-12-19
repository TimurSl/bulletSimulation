import requests
import math


def get_wind_vector(api_key, location=None, coordinates=None):
    """
    Get wind speed in all axes (X, Y) based on location or coordinates.

    Parameters:
        api_key (str): Your OpenWeatherMap API key.
        location (str): City name (e.g., "London").
        coordinates (tuple): Latitude and longitude as (lat, lon).

    Returns:
        dict: Wind speed in m/s for X, Y and additional metadata.
    """
    # Base URL for OpenWeatherMap
    base_url = "https://api.openweathermap.org/data/2.5/weather"

    # Prepare query parameters
    params = {"appid": api_key, "units": "metric"}
    if location:
        params["q"] = location
    elif coordinates:
        params["lat"] = coordinates[0]
        params["lon"] = coordinates[1]
    else:
        raise ValueError("Either 'location' or 'coordinates' must be provided.")

    # Send request to OpenWeatherMap
    response = requests.get(base_url, params=params)
    if response.status_code != 200:
        raise Exception(f"API request failed: {response.json().get('message', 'Unknown error')}")

    # Extract wind data
    data = response.json()
    wind_speed = data["wind"]["speed"]  # Wind speed in m/s
    wind_deg = data["wind"]["deg"]  # Wind direction in degrees (meteorological)

    # Convert meteorological direction to cartesian vector (X, Y)
    wind_rad = math.radians(wind_deg)
    wind_x = wind_speed * math.sin(wind_rad)  # Wind component along X-axis
    wind_y = wind_speed * math.cos(wind_rad)  # Wind component along Y-axis

    return {
        "wind_x": wind_x,
        "wind_y": wind_y,
        "wind_speed": wind_speed,
        "wind_direction": wind_deg,
        "location": data["name"],
        "coordinates": (data["coord"]["lat"], data["coord"]["lon"]),
    }


# Example usage
if __name__ == "__main__":
    api_key = "c3c167faed44acf64accc2ac70b1b3c5"  # Replace with your OpenWeatherMap API key

    # Example: Query by city name
    location = "Kyiv"
    result = get_wind_vector(api_key, location=location)
    print(f"Wind in {result['location']} (Coordinates: {result['coordinates']}):")
    print(f"  - Speed: {result['wind_speed']} m/s")
    print(f"  - Direction: {result['wind_direction']}°")
    print(f"  - X-axis component: {result['wind_x']} m/s")
    print(f"  - Y-axis component: {result['wind_y']} m/s")

    # Example: Query by coordinates
    coordinates = (51.5074, -0.1278)  # Latitude, Longitude for London
    result = get_wind_vector(api_key, coordinates=coordinates)
    print(f"\nWind at coordinates {result['coordinates']}:")
    print(f"  - Speed: {result['wind_speed']} m/s")
    print(f"  - Direction: {result['wind_direction']}°")
    print(f"  - X-axis component: {result['wind_x']} m/s")
    print(f"  - Y-axis component: {result['wind_y']} m/s")
