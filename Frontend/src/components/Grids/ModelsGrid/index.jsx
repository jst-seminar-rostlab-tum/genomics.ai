import React from 'react';
import { Box } from '@mui/material';
import Grid from '@mui/material/Grid';
import { ModelCard } from 'components/Cards/ModelCard';

import styles from './modelsGrid.module.css';

const ModelsGrid = ({ models, path, compatibleModels = null, selectedModel = null, handleModelSelection = null }) => {

  console.log(compatibleModels)

  return (
    <Box className={styles.cardsContainer} mb="2em" >
      <Grid container spacing={3}>
        {models && models.map((model) => (
          <Grid key={model._id} item xs={12} sm={6} md={4} lg={3}>
            <ModelCard
              title={`Model ${model.name}`}
              description={model.description}
              learnMoreLink={`${path}/models/${model._id}`}
              disabled={compatibleModels && (!compatibleModels.map((m) => m.toLowerCase()).includes(model.name.toLowerCase()) || compatibleModels.length == 0)}
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
