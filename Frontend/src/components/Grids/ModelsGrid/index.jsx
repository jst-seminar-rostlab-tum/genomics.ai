import React from 'react';
import { Box } from '@mui/material';
import Grid from '@mui/material/Grid';
import { ModelCard } from 'components/Cards/ModelCard';
import { useState } from 'react';

import styles from './modelsGrid.module.css';

const ModelsGrid = ({ models, path, selectedAtlas = null, selectedModel = null, handleModelSelection = null }) => {

  const checkIfDisabled = (name) => {
    return !compatibleModels.map(m => m.toLowerCase()).includes(m.name.toLowerCase()) || compatibleModels.length == 0
  }

  return (
    <Box className={styles.cardsContainer} maxHeight="50vh" mb="2em" >
      <Grid container spacing={3}>
        {models && models.map((model) => (
          <Grid key={model._id} item xs={12} sm={6} md={4} lg={3}>
            <ModelCard
              title={`Model ${model.name}`}
              description={model.description}
              learnMoreLink={`${path}/models/${model._id}`}
              disabled={!compatibleModels.map(m => m.toLowerCase()).includes(m.name.toLowerCase()) || compatibleModels.length == 0}
              onSelect={() => { if(handleModelSelection) handleModelSelection(model) }}
              selected={selectedModel && selectedModel.name === model.name}
            />
          </Grid>
        ))}
      </Grid>
    </Box>
  )
}

export default ModelsGrid;
