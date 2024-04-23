from django.core.management.base import BaseCommand
from vetu.models import Author  # Import your Author model
import openai

class Command(BaseCommand):
    def handle(self, *args, **options):
        # Set your OpenAI API key
        openai.api_key = 'YOUR_API_KEY'

        # Fetch data from the Author model
        authors = Author.objects.all()  # Fetch all authors or apply your filtering logic

        # Prepare data for prompt
        for author in authors:
            affiliation = author.affiliations
            decoding_dictionary = {
                # Define your decoding dictionary based on your requirements
                # Example: '1': 'Surgical Ward', '2': 'Medical Ward', '3': 'Intensive Care Unit', etc.
            }

            # Format the prompt for each affiliation
            prompt = f"For the following affiliation '{affiliation}', which of the following hospital ward types is the most apt? {decoding_dictionary}\nAI:"

            # Call the OpenAI API to get a response
            response = openai.Completion.create(
                engine="text-davinci-003-turbo",  # Using ChatGPT 3.5 Turbo
                prompt=prompt,
                max_tokens=150
            )

            # Extract and display the AI's reply
            ai_reply = response['choices'][0]['text']
            self.stdout.write(self.style.SUCCESS(f"AI Response for affiliation '{affiliation}': {ai_reply}"))
