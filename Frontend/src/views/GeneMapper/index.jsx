import { Box, Button } from '@mui/material';
import { useCallback } from 'react';
import { Redirect, Route, Switch, useHistory, useRouteMatch } from 'react-router-dom';
import GeneMapperState from './GeneMapperState';
import GeneMapperHome from './Home';
import GeneMapperResultView from './Result';

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
          <Button onClick={() => history.push(`${path}/create`)}>Create</Button>
          <Button onClick={() => history.push(`${path}/result`)}>Result</Button>
          <GeneMapperHome basePath={path} />
        </Route>

        <Route path={`${path}/create`}>
          <GeneMapperState basePath={`${path}/create`} path={`${path}`} />
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
