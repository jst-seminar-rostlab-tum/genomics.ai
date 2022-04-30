import React, { useCallback } from 'react';
import { Box, Button } from '@mui/material';
import {
  Switch, Route, Redirect, useHistory, useRouteMatch,
} from 'react-router-dom';
import GeneMapperResultView from './Result';
import AtlasModelChoice from './AtlasModelChoice/AtlasModelChoice';
import UploadFilePage from './UploadFilePage';
import GeneMapperHome from './Home';

function GeneMapper({ sidebarShown }) {
  const paddingL = useCallback(() => (sidebarShown ? '100px' : '350px'), [sidebarShown]);

  const { path, url } = useRouteMatch();
  const history = useHistory();

  return (
    <Box sx={{
      pl: paddingL,
      pr: '20px',
      height: '100vh',
      width: '100vw',
    }}
    >
      <Switch>
        <Route exact path={`${path}`}>
          <Button onClick={() => history.push(`${path}/selection`)}>Selection</Button>
          <Button onClick={() => history.push(`${path}/upload`)}>Upload</Button>
          <Button onClick={() => history.push(`${path}/result`)}>Result</Button>
          <GeneMapperHome basePath={path} />
        </Route>

        <Route path={`${path}/selection`}>
          <AtlasModelChoice />
        </Route>

        <Route path={`${path}/upload`}>
          <UploadFilePage basePath={`${path}/upload`} />
        </Route>

        <Route path={`${path}/result/:projectId`}>
          <GeneMapperResultView />
        </Route>

        <Route path={`${path}/*`}>
          <Redirect to={`${url}`} />
        </Route>
      </Switch>
    </Box>
  );
}

export default GeneMapper;
