# Choose Atlas and Model - View

This view shows the first page of a two step guide for creating a new mapping.


## Main functionalities of this view:

* The user selects the atlas to be used for the mapping
* Next the user selects the model to be used for the mapping
* By clicking the Cancel-Button on the bottom left the user gets redirected to the GeneMapper-Homeview
* By clicking the Confirm-Button on the bottom right the user gets forwarded to page 2 of the 2-step guide


## Additional information:

The models that can be used for the mapping depend on the selected atlas. That is why when no atlas is selected all models are greyed out and disabled. After selecting an atlas the compatible models get selectable. 
You can only click confirm when an atlas and a model where selected. When clicking confirm even though a selection is missing, a warning shows up.


## About the code

* The shown atlases and models as well as a boolean whether the Confirm-Warning is active is stored as a local state in the AtlasModelChoice functional component.
* The useEffect hook is used to fetch the atlases and models at first render while also adjusting their data for rendering
