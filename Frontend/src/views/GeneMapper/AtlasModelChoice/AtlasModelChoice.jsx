import {
  Grid, Typography, Stack, Alert, Box, Tooltip,
} from '@mui/material';
import AtlasCardSelect from 'components/Cards/AtlasCardSelect';
import { TabCard } from 'components/GeneMapper/TabCard';
import { ModelCardSelect } from 'components/Cards/ModelCardSelect';
import CustomButton from 'components/CustomButton';
import { colors } from 'shared/theme/colors';
import { useHistory } from 'react-router-dom';
import React, { useState, useEffect } from 'react';
import Clear from '@mui/icons-material/Clear';
import ArrowForwardIcon from '@mui/icons-material/ArrowForward';

function AtlasModelChoice({
  setActiveStep,
  selectedAtlas, setSelectedAtlas,
  selectedModel, setSelectedModel, path,
  compatibleModels, atlases, models,
  selectedDataset, setSelectedDataset,
  datasetIsSelected, setDatasetIsSelected,
  demos,
}) {
  const [showWarning, setShowWarning] = useState(false);
  const history = useHistory();

  /* Handling the choice of the demo dataset */
  const handleDemoClick = (dataset) => {
    // only set selected dataset if not already selected or null
    if (!selectedDataset || selectedDataset._id !== dataset._id) {
      setSelectedDataset(dataset);
      setDatasetIsSelected(true);
      // find the objects corresponding to the atlas and model in the array
      const atlasObj = atlases.filter(
        (a) => a.name.toLowerCase() === dataset.atlas.toLowerCase(),
      )[0];
      const modelObj = models.filter(
        (m) => m.name.toLowerCase() === dataset.model.toLowerCase(),
      )[0];

      setSelectedAtlas(atlasObj);
      setSelectedModel(modelObj);
    } else {
      setDatasetIsSelected(false);
      setSelectedDataset(null);
      // deselect the atlas and model
      setSelectedAtlas({});
      setSelectedModel({});
    }
  };

  return (
    <div>
      {showWarning
        && (
          <Alert severity="error">
            Select an Atlas and a fitting Model before continuing
          </Alert>
        )}
      <Typography
        variant="h5"
        sx={{
          fontWeight: 'bold',
          pb: '1em',
          mt: '1.5em',
        }}
      >
        Pick an Atlas
      </Typography>

      <Grid container spacing={2} width="100%" overflow="auto" wrap="nowrap">
        {console.log(atlases)}
        {
          atlases && atlases.map((a) => (
            <Grid item height="320px">
              <AtlasCardSelect
                width="225px"
                height="97%"
                title={a.name}
                modalities={a.modalities}
                cellsInReference={a.numberOfCells}
                species={a.species}
                imgLink={a.previewPictureURL}
                selected={selectedAtlas.name === a.name}
                onSelect={setSelectedAtlas}
                atlasObject={a}
              />
            </Grid>
          ))
        }
      </Grid>

      <Typography variant="h5" sx={{ fontWeight: 'bold', pb: '1em' }} marginTop="32px">
        Pick a Model
      </Typography>

      <Grid container spacing={2} width="100%" overflow="auto" wrap="nowrap">
        {
          models && models.map((m) => (
            <Grid item height="225px">
              <ModelCardSelect
                width="225px"
                height="97%"
                title={m.name}
                description={m.description}
                selected={selectedModel.name === m.name}
                onSelect={setSelectedModel}
                modelObject={m}
                disabled={!compatibleModels
                  || !compatibleModels.map(
                    (m) => m.toLowerCase(),
                  ).includes(m.name.toLowerCase()) || compatibleModels.length === 0}
              />
            </Grid>
          ))
        }
      </Grid>

      {/* Demo datasets */}
      <Box sx={{ width: 400 }}>
        <Typography variant="h5" sx={{ fontSize: '18px', fontWeight: 'bold', pb: '1em' }} marginTop="32px">
          Or try out one of the available demo datasets
        </Typography>
      </Box>
      {/* Loop over demo datasets */}
      { demos && demos.map((dataset) => (
        <TabCard
          width="40%"
          height="50px"
          title={dataset.name}
          atlas={dataset.atlas}
          model={dataset.model}
          handleOnClick={() => handleDemoClick(dataset)}
          selected={datasetIsSelected && dataset._id === selectedDataset._id}
        />
      ))}

      <Stack direction="row" justifyContent="space-between" sx={{ marginTop: '50px', marginBottom: '3em' }}>
        <CustomButton type="tertiary" onClick={() => history.push(`${path}`)}>
          <Clear />
          &nbsp; Cancel
        </CustomButton>
        <Tooltip title={!selectedAtlas || !selectedModel ? 'Select an Atlas and a fitting Model before continuing' : ''} placement="top">
          <Box
            onClick={!selectedAtlas || !selectedModel ? () => setShowWarning(true) : () => { }}
          >
            <CustomButton
              type="primary"
              disabled={!selectedAtlas || !selectedModel}
              onClick={() => setActiveStep(1)}
            >
              Confirm&nbsp;
              <ArrowForwardIcon />
            </CustomButton>
          </Box>
        </Tooltip>
      </Stack>
    </div>
  );
}

export default AtlasModelChoice;
