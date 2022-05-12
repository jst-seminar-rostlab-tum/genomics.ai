import React from 'react';
import { Box } from '@mui/material';
import Grid from '@mui/material/Grid';
import { ModelCard } from 'components/Cards/ModelCard';
import { useLocation } from 'react-router-dom';

import styles from './modelsGrid.module.css';

function applyModelFilters(models, searchedKeyword, searchParams) {
  const searchedModels = models.filter(
    (item) => item.name.toLowerCase().includes(searchedKeyword.toLowerCase()),
  );
  if (searchParams.get('sortBy') === 'name' || searchParams.get('sortBy') === null) {
    searchedModels.sort((a, b) => {
      const nameA = a.name.toUpperCase();
      const nameB = b.name.toUpperCase();
      if (nameA < nameB) {
        return -1;
      }
      if (nameA > nameB) {
        return 1;
      }
      return 0;
    });
  }
  return searchedModels;
}

const ModelsGrid = ({ models, searchedKeyword, path }) => {
  const { search } = useLocation();
  const searchParams = new URLSearchParams(search);
  const modelsFiltered = applyModelFilters(models, searchedKeyword, searchParams);
  return (
    <Box className={styles.cardsContainer}>
      <Grid container spacing={3}>
        {modelsFiltered && modelsFiltered.map((model) => (
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
};
export default ModelsGrid;
