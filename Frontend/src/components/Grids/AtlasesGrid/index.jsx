import React from 'react';
import { useHistory } from 'react-router-dom';
import { Box } from '@mui/material';
import Grid from '@mui/material/Grid';
import AtlasCard from 'components/Cards/AtlasCard';
import styles from './atlasesGrid.module.css';

const AtlasesGrid = ({ atlases, path }) => {
  const history = useHistory();
  return (
    <Box className={styles.atlasContainer}>
      <Grid container spacing={3}>
        {atlases && atlases.map((atlas) => (
          <Grid key={atlas._id} item xs={12} sm={6} md={4} lg={3}>
            <AtlasCard
              onClick={() => history.push(`${path}/atlases/${atlas._id}/visualization`)}
              atlasId={atlas._id}
              imgLink={atlas.previewPictureURL}
              species={atlas.species}
              modalities={atlas.modalities}
              title={atlas.name}
              learnMoreLink={`${path}/atlases/${atlas._id}`}
            />
          </Grid>
        ))}
        {atlases && atlases.map((atlas) => (
          <Grid key={atlas._id} item xs={12} sm={6} md={4} lg={3}>
            <AtlasCard
              onClick={() => history.push(`${path}/atlases/${atlas._id}/visualization`)} // removed setSelectedAtlas => feels outdated
              atlasId={atlas._id}
              imgLink={atlas.previewPictureURL}
              species={atlas.species}
              modalities={atlas.modalities}
              title={atlas.name}
              learnMoreLink={`${path}/atlases/${atlas._id}`}
            />
          </Grid>
        ))}
        {atlases && atlases.map((atlas) => (
          <Grid key={atlas._id} item xs={12} sm={6} md={4} lg={3}>
            <AtlasCard
              onClick={() => history.push(`${path}/atlases/${atlas._id}/visualization`)} // removed setSelectedAtlas => feels outdated
              atlasId={atlas._id}
              imgLink={atlas.previewPictureURL}
              species={atlas.species}
              modalities={atlas.modalities}
              title={atlas.name}
              learnMoreLink={`${path}/atlases/${atlas._id}`}
            />
          </Grid>
        ))}
      </Grid>
    </Box>
  );
};

export default AtlasesGrid;
