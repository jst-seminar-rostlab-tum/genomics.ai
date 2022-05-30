# Upload file - View

This view show the second page of a two step guide for creating a new mapping.


## Main functionalities of this view:

* It shows a summary of the selections taken at the previous step
* The consequent requirements that result from the selections are displayed. These are relevant information that need to be considered for the file upload.
* The user is able to select a file for the mapping by browsing or using drag-and-drop. Or existing datasets from previous mapping can be used.
* By clicking the Back-Button on the bottom left the user gets redirected to the previous step.
* By clicking the Create-Mapping-Button the user gets forwarded to a final popup where he can name and submit the mapping-project.
* By clicking Learn-More on the Atlas or Model Card a popup with further information shows up.


## About the code:

* In the useEffect hook the requirements that depend on the chosen model are fetched and stored in a local state
* local states and helper functions are described in detail in the comments in the code.