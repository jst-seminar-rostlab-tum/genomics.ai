import React from 'react';
import { Box } from '@mui/material';
import Grid from '@mui/material/Grid';
import { ModelCard } from 'components/Cards/ModelCard';

import styles from './modelsGrid.module.css';

const ModelsGrid = ({ models, path }) => (
  <Box className={styles.cardsContainer} maxHeight="50vh" mb="2em">
    <Grid container spacing={3}>
      {models && models.map((model) => (
        <Grid key={model._id} item xs={12} sm={6} md={4} lg={3}>
          <ModelCard
            onClick={() => {}} // removed setSelectedModel => feels outdated
            title={`Model ${model.name}`}
            description={model.description}
            learnMoreLink={`${path}/models/${model._id}`}
          />
        </Grid>
      ))}
    </Grid>
  </Box>
);
export default ModelsGrid;
